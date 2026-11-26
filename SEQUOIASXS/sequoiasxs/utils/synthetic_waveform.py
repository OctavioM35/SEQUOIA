import numpy as np
import bilby
import sxs
from bilby.gw.detector import InterferometerList
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from gwsurrogate.surrogate import _sYlm
from sequoiasxs.waveforms.dansur.waveform_generator import sur, _numba_taper_modes, waveform_generator


def _get_float(d, key, default):
    val = d.get(key, default)
    try:
        return float(val)
    except (ValueError, TypeError):
        return float(default)


def synthetic_waveform(dicc):
    """
    Generate synthetic gravitational-wave detector data (H1 and L1) from an
    SXS numerical-relativity waveform, strictly aligned with DANSur's cyclic
    time convention, domain, taper, modal phases, and Bilby's frequency-domain projection.
    """
    # ---------------------------------------------------------
    # 1. Configuration & defaults from dicc
    # ---------------------------------------------------------
    total_mass = _get_float(dicc, "total_mass", 70.0)                         # M_sun
    luminosity_distance_mpc = _get_float(dicc, "luminosity_distance", 50.0)  # Mpc
    fs = int(_get_float(dicc, "sampling-frequency", 4096))                    # Hz
    psd_value = _get_float(dicc, "psd_value", 1e-45)                          # flat PSD [1/Hz]

    # Angles and sky position
    theta = _get_float(dicc, "theta", 0.9)      # inclination / polar angle (theta_jn)
    phi = _get_float(dicc, "phi", 0.0)          # coalescence phase (Bilby convention)
    ra = _get_float(dicc, "ra", 1.375)          # right ascension
    dec = _get_float(dicc, "dec", -1.21)        # declination
    psi = _get_float(dicc, "psi", 0.0)          # polarization angle

    default_modes = [(2, 2), (3, 3), (2, 1), (4, 4), (3, 2), (5, 5), (4, 3)]
    modes = dicc.get("modes", default_modes)

    # Physical constants (SI units)
    G = 6.67430e-11
    M_sun = 1.98847e30
    c = 299792458.0
    Mpc = 3.085677581e22
    M_sun_time = 4.9254909476412675e-6
    time_conversion = total_mass * M_sun_time

    # ---------------------------------------------------------
    # 2. Load SXS simulation & extract metadata
    # ---------------------------------------------------------
    sxs_id = dicc["sxs"]
    try:
        sim_obj = sxs.load(sxs_id, ignore_deprecation=True)
    except Exception as e:
        print(e)
        sim_obj = sxs.load(sxs_id)

    h = sim_obj.h if hasattr(sim_obj, "h") else sim_obj

    m1_ref = float(h.metadata.get("reference_mass1"))
    m2_ref = float(h.metadata.get("reference_mass2"))
    if m1_ref is None or m2_ref is None:
        raise ValueError(f"Missing reference masses in SXS metadata for {sxs_id}.")

    spin1_ref = h.metadata.get("reference_dimensionless_spin1", [0.0, 0.0, 0.0])
    spin2_ref = h.metadata.get("reference_dimensionless_spin2", [0.0, 0.0, 0.0])

    # Bloqueo de precesión real (umbral físico 1e-3) y órbitas no circulares
    in_plane_1 = float(np.hypot(spin1_ref[0], spin1_ref[1])) if hasattr(spin1_ref, "__getitem__") else 0.0
    in_plane_2 = float(np.hypot(spin2_ref[0], spin2_ref[1])) if hasattr(spin2_ref, "__getitem__") else 0.0
    if max(in_plane_1, in_plane_2) >= 1e-3:
        raise ValueError(
            f"{sxs_id} has precession (spin_xy = {max(in_plane_1, in_plane_2):.2e} >= 1e-3). Rejected"
        )

    ecc_raw = h.metadata.get("reference_eccentricity", None)
    ecc_str = str(ecc_raw).strip().lstrip("<") if ecc_raw is not None else "nan"
    try:
        ecc_val = float(ecc_str)
    except ValueError:
        ecc_val = np.inf
    if not np.isfinite(ecc_val) or ecc_val > 2e-3:
        raise ValueError(
            f"{sxs_id} is not quasi-circular (reference_eccentricity = {ecc_raw}). Rejected."
        )

    # Ensure convention m1 >= m2 (q = m2/m1 <= 1) and swap spins accordingly
    if m1_ref >= m2_ref:
        mass1 = total_mass * m1_ref / (m1_ref + m2_ref)
        mass2 = total_mass * m2_ref / (m1_ref + m2_ref)
        spin1, spin2 = spin1_ref, spin2_ref
    else:
        mass1 = total_mass * m2_ref / (m1_ref + m2_ref)
        mass2 = total_mass * m1_ref / (m1_ref + m2_ref)
        spin1, spin2 = spin2_ref, spin1_ref

    q_inj = float(mass2 / mass1)
    chirp_mass = float((mass1 * mass2) ** 0.6 / (mass1 + mass2) ** 0.2)

    chi_1 = float(spin1[2] if hasattr(spin1, "__getitem__") else spin1)
    chi_2 = float(spin2[2] if hasattr(spin2, "__getitem__") else spin2)

    # Reference orbital frequency
    f_orb_init = h.metadata.get("initial_orbital_frequency", None)
    if f_orb_init is None:
        f_orb_init = h.metadata.get("reference_orbital_frequency", 0.015)
    f_orb_val = float(f_orb_init[2] if hasattr(f_orb_init, "__getitem__") else f_orb_init)

    ref_orb_freq = h.metadata.get("reference_orbital_frequency", f_orb_val)
    ref_orb_z = float(ref_orb_freq[2] if hasattr(ref_orb_freq, "__getitem__") else ref_orb_freq)
    f_ref = ref_orb_z * c**3 / (np.pi * G * total_mass * M_sun)

    # ---------------------------------------------------------
    # 3. Filter modes, match DANSur domain & trim valid SXS window
    # ---------------------------------------------------------
    dansur_query = {
        "q": float(q_inj),
        "chiA0": float(chi_1),
        "chiB0": float(chi_2),
        "M": float(total_mass),
        "dist_mpc": float(luminosity_distance_mpc),
        "f_low": 0,
        "times": None,
        "f_ref": None,
        "freqs": None,
        "inclination": None,
        "phi_ref": 0.0,
        "units": "mks",
    }
    domain_dansur, h_modes_dan, _ = sur(**dansur_query)
    t_start_dansur_s = float(domain_dansur[0])
    t_end_dansur_s = float(domain_dansur[-1])

    # Filter exact (l, +-m) modes in SXS
    target_modes = set()
    for l_val, m_val in modes:
        target_modes.add((int(l_val), int(m_val)))
        target_modes.add((int(l_val), -int(m_val)))

    keep_indices = [
        i for i, (l_val, m_val) in enumerate(h.LM)
        if (int(l_val), int(m_val)) in target_modes
    ]
    if not keep_indices:
        raise ValueError(f"None of the requested modes {modes} were found in {sxs_id}.")

    h_data = np.array(h.ndarray, copy=True)
    keep_mask = np.zeros(len(h.LM), dtype=bool)
    keep_mask[keep_indices] = True
    h_data[:, ~keep_mask] = 0.0

    # Align t=0 with the L2 norm across active modes
    t_all = np.asarray(h.t)
    mode_norm = np.linalg.norm(h_data, axis=1)
    t_peak = t_all[np.argmax(mode_norm)]

    relaxation = float(h.metadata.get("relaxation_time", t_all[0]))
    t_start_sxs_s = (relaxation - t_peak) * time_conversion

    t_inspiral_start_s = max(t_start_sxs_s, t_start_dansur_s)
    t_ringdown_end_s = min((t_all[-1] - t_peak) * time_conversion, t_end_dansur_s)

    t_seconds_all = (t_all - t_peak) * time_conversion
    valid_sxs = (t_seconds_all >= t_inspiral_start_s) & (t_seconds_all <= t_ringdown_end_s)

    # ---------------------------------------------------------
    # 4. Dynamic buffer timing & duration
    # ---------------------------------------------------------
    t_inspiral = abs(t_inspiral_start_s)
    t_ringdown = max(0.1, t_ringdown_end_s)

    min_required_duration = t_inspiral + t_ringdown + 0.6
    if "duration" in dicc and isinstance(dicc["duration"], (int, float)) and dicc["duration"] >= min_required_duration:
        duration = float(dicc["duration"])
    else:
        duration = float(max(2.0, 2 ** int(np.ceil(np.log2(min_required_duration)))))

    N = int(duration * fs)
    N -= N % 2

    t_merger_offset = duration - max(0.3, t_ringdown + 0.15)
    gt_val = dicc.get("geocent_time", 1126259462.0)
    geocent_time = float(gt_val) if isinstance(gt_val, (int, float)) else 1126259462.0
    start_time = geocent_time - t_merger_offset

    dicc["geocent_time"] = float(geocent_time)
    dicc["start_time"] = float(start_time)
    dicc["duration"] = float(duration)
    dicc["sampling-frequency"] = int(fs)
    dicc["t_merger_offset"] = float(t_merger_offset)

    # ---------------------------------------------------------
    # 5. Interpolate onto uniform linear grid, taper & auto-align with DANSur
    # ---------------------------------------------------------
    t_src = t_seconds_all[valid_sxs]
    t_uniform = np.arange(t_src[0], t_src[-1], 1.0 / fs)
    n_u = len(t_uniform)
    n_active = len(keep_indices)

    modes_u_r = np.zeros((n_active, n_u), dtype=np.float64)
    modes_u_i = np.zeros((n_active, n_u), dtype=np.float64)

    for k, idx_mode in enumerate(keep_indices):
        hlm_src = h_data[valid_sxs, idx_mode]
        modes_u_r[k] = np.interp(t_uniform, t_src, hlm_src.real)
        modes_u_i[k] = np.interp(t_uniform, t_src, hlm_src.imag)

    # Apply Planck taper at the chronological start of the inspiral
    modes_u_r, modes_u_i = _numba_taper_modes(modes_u_r, modes_u_i)

    # Estimate f_min right after taper completion on uniform grid
    idx_22_local = [k for k, idx_m in enumerate(keep_indices) if tuple(h.LM[idx_m]) == (2, 2)][0]
    h22_u = modes_u_r[idx_22_local] + 1j * modes_u_i[idx_22_local]
    amp_22 = np.abs(h22_u)
    taper_end = np.where(amp_22 >= 0.95 * np.max(amp_22[:max(10, n_u // 4)]))[0]
    idx_fstart = int(taper_end[0]) + 5 if len(taper_end) > 0 else int(0.1 * n_u)
    idx_fstart = min(idx_fstart, n_u - 2)

    phase_22 = np.unwrap(np.angle(h22_u))
    dphi_dt = np.gradient(phase_22, t_uniform)
    f_start_trimmed = abs(float(dphi_dt[idx_fstart])) / (2.0 * np.pi)
    f_min = float(max(20.0, f_start_trimmed * 1.25))

    dicc.setdefault("waveform", {})
    dicc["waveform"]["minimum-frequency"] = float(f_min)
    dicc["waveform"]["reference-frequency"] = float(f_ref)
    dicc["waveform"]["f_low"] = float(f_min)
    dicc["waveform"]["f_ref"] = float(f_ref)

    factor = (G * total_mass * M_sun) / (c**2 * luminosity_distance_mpc * Mpc)

    # 5.1 Multi-branch & parity alignment directly maximizing frequency-domain overlap
    t_align_sec = -100.0 * time_conversion
    idx_align_sxs = int(np.argmin(np.abs(t_uniform - t_align_sec)))
    idx_align_dan = int(np.argmin(np.abs(np.asarray(domain_dansur) - t_align_sec)))

    h22_sxs_val = modes_u_r[idx_22_local, idx_align_sxs] + 1j * modes_u_i[idx_22_local, idx_align_sxs]
    h22_dan_val = h_modes_dan[(2, 2)][idx_align_dan]
    base_dphi = 0.5 * float(np.angle(h22_sxs_val * np.conj(h22_dan_val)))

    def _synthesize_and_fft(delta_phi_orb=0.0, delta_t=0.0, global_sign=1.0, parity_sign=1):
        hp_u = np.zeros(n_u, dtype=np.float64)
        hc_u = np.zeros(n_u, dtype=np.float64)
        for k, idx_mode in enumerate(keep_indices):
            l_val, m_val = h.LM[idx_mode]
            parity = float(parity_sign ** abs(int(m_val)))
            rot = global_sign * parity * np.exp(-1j * int(m_val) * delta_phi_orb)
            hlm_rot = (modes_u_r[k] + 1j * modes_u_i[k]) * rot
            ylm = _sYlm(-2, int(l_val), int(m_val), float(theta), float(-phi))
            hr = ylm.real * hlm_rot.real - ylm.imag * hlm_rot.imag
            hi = ylm.real * hlm_rot.imag + ylm.imag * hlm_rot.real
            hp_u += factor * hr
            hc_u += -factor * hi

        t_elapsed = np.arange(N) / fs
        t_post_max = t_ringdown + 0.05
        tau_cyclic = np.where(t_elapsed <= t_post_max, t_elapsed, t_elapsed - duration) - delta_t

        hp_cyc = np.interp(tau_cyclic, t_uniform, hp_u, left=0.0, right=0.0)
        hc_cyc = np.interp(tau_cyclic, t_uniform, hc_u, left=0.0, right=0.0)
        hp_f, _ = bilby.core.utils.nfft(hp_cyc, fs)
        hc_f, _ = bilby.core.utils.nfft(hc_cyc, fs)
        return hp_u, hc_u, {"plus": hp_f, "cross": hc_f}

    injection_parameters = {
        "mass_1": float(mass1),
        "mass_2": float(mass2),
        "mass_ratio": float(q_inj),
        "chirp_mass": float(chirp_mass),
        "luminosity_distance": float(luminosity_distance_mpc),
        "chi_1": float(chi_1),
        "chi_2": float(chi_2),
        "theta_jn": float(theta),
        "phase": float(phi),
        "ra": float(ra),
        "dec": float(dec),
        "psi": float(psi),
        "geocent_time": float(geocent_time),
    }

    freqs = np.linspace(0, fs / 2, N // 2 + 1)
    psd_array = np.ones_like(freqs) * psd_value
    ifo_ref = bilby.gw.detector.get_empty_interferometer("H1")
    ifo_ref.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
        frequency_array=freqs, psd_array=psd_array
    )
    ifo_ref.minimum_frequency = float(f_min)
    ifo_ref.maximum_frequency = float(fs / 2)
    ifo_ref.strain_data.roll_off = 0.0
    ifo_ref.set_strain_data_from_zero_noise(
        sampling_frequency=fs, duration=duration, start_time=start_time
    )

    wf_gen = waveform_generator(dicc, list(injection_parameters.keys()))
    pols_dansur = wf_gen.frequency_domain_strain(injection_parameters)
    h_dansur_fd = ifo_ref.get_detector_response(pols_dansur, injection_parameters)

    mask_f = (ifo_ref.frequency_array >= f_min) & (ifo_ref.frequency_array <= fs / 2)
    psd_m = ifo_ref.power_spectral_density_array[mask_f]
    df = ifo_ref.frequency_array[1] - ifo_ref.frequency_array[0]
    N_full = (len(h_dansur_fd) - 1) * 2

    # 5.2 Search over all 8 (parity, branch) combinations using full waveform correlation
    best_max_corr = -1.0
    best_dphi = base_dphi
    best_parity_sign = 1
    best_dt_align = 0.0
    best_global_sign = 1.0

    for p_sign in (1, -1):
        for k_branch in range(4):
            cand_dphi = base_dphi + k_branch * (np.pi / 2.0)
            _, _, pols_cand = _synthesize_and_fft(
                delta_phi_orb=cand_dphi, delta_t=0.0, global_sign=1.0, parity_sign=p_sign
            )
            h_cand_fd = ifo_ref.get_detector_response(pols_cand, injection_parameters)

            integrand = 4.0 * df * (h_cand_fd[mask_f] * np.conj(h_dansur_fd[mask_f]) / psd_m)
            full_pos = np.zeros(N_full, dtype=np.complex128)
            full_pos[:len(mask_f)][mask_f] = integrand
            corr = np.fft.ifft(full_pos) * N_full
            abs_corr = np.abs(corr)
            idx_b = int(np.argmax(abs_corr))

            if abs_corr[idx_b] > best_max_corr:
                best_max_corr = float(abs_corr[idx_b])
                best_dphi = cand_dphi
                best_parity_sign = p_sign

                ym1 = float(abs_corr[(idx_b - 1) % N_full])
                y0 = float(abs_corr[idx_b])
                yp1 = float(abs_corr[(idx_b + 1) % N_full])
                denom = ym1 - 2.0 * y0 + yp1
                delta_idx = float(np.clip(0.5 * (ym1 - yp1) / denom, -0.5, 0.5)) if abs(denom) > 1e-30 else 0.0
                idx_signed = (idx_b if idx_b < N_full // 2 else idx_b - N_full) + delta_idx
                best_dt_align = idx_signed / fs
                best_global_sign = -1.0 if np.real(corr[idx_b]) < 0 else 1.0

    # 5.2b Continuous refinement of orbital phase so Direct Overlap matches Max Overlap
    def _neg_direct_overlap(dphi_offset):
        _, _, pols_try = _synthesize_and_fft(
            delta_phi_orb=best_dphi + dphi_offset,
            delta_t=-best_dt_align,
            global_sign=best_global_sign,
            parity_sign=best_parity_sign,
        )
        h_try = ifo_ref.get_detector_response(pols_try, injection_parameters)
        num = 4.0 * df * np.real(np.sum(h_try[mask_f] * np.conj(h_dansur_fd[mask_f]) / psd_m))
        norm_try = np.sqrt(4.0 * df * np.real(np.sum(np.abs(h_try[mask_f])**2 / psd_m)))
        norm_dan = np.sqrt(4.0 * df * np.real(np.sum(np.abs(h_dansur_fd[mask_f])**2 / psd_m)))
        return -float(num / (norm_try * norm_dan))

    res_phi = minimize_scalar(_neg_direct_overlap, bounds=(-np.pi / 4.0, np.pi / 4.0), method="bounded")
    best_dphi = float(best_dphi + res_phi.x)

    # 5.3 Final exact synthesis with verified alignment parameters
    h_plus_sxs, h_cross_sxs, polarizations = _synthesize_and_fft(
        delta_phi_orb=best_dphi,
        delta_t=-best_dt_align,
        global_sign=best_global_sign,
        parity_sign=best_parity_sign,
    )
    dt_align = best_dt_align
    t_seconds = t_uniform

    # ---------------------------------------------------------
    # 6. Initialize interferometers & inject in frequency domain
    # ---------------------------------------------------------
    interferometers = []
    signals = {}

    for name in ("H1", "L1"):
        ifo = bilby.gw.detector.get_empty_interferometer(name)
        ifo.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
            frequency_array=freqs, psd_array=psd_array
        )
        ifo.minimum_frequency = float(f_min)
        ifo.maximum_frequency = float(fs / 2)

        # Disable Tukey roll_off to match WaveformGenerator's raw FFT
        ifo.strain_data.roll_off = 0.0
        ifo.set_strain_data_from_zero_noise(
            sampling_frequency=fs, duration=duration, start_time=start_time
        )
        ifo.strain_data.roll_off = 0.0

        ifo.inject_signal(
            parameters=injection_parameters,
            injection_polarizations=polarizations,
        )

        sig_fd = ifo.get_detector_response(polarizations, injection_parameters)
        signals[name] = sig_fd
        interferometers.append(ifo)

    ifos = InterferometerList(interferometers)

    # ---------------------------------------------------------
    # 7. Store injection dictionary & run validation diagnostics
    # ---------------------------------------------------------
    signal_h1_td = bilby.core.utils.infft(signals["H1"], fs)
    dicc["signal"] = [signal_h1_td, h_cross_sxs, h_plus_sxs]
    dicc["inyection"] = injection_parameters
    dicc["chirp_mass"] = float(chirp_mass)

    print("=== Injection Summary & Validation Tests ===")
    snr_dict = {}
    for ifo in ifos:
        snr_sq = ifo.optimal_snr_squared(signals[ifo.name])
        snr_dict[ifo.name] = float(np.sqrt(np.real(snr_sq)))
        print(f"  {ifo.name} optimal SNR          : {snr_dict[ifo.name]:.4f}")
    snr_net = np.sqrt(sum(s**2 for s in snr_dict.values()))
    print(f"  Network optimal SNR     : {snr_net:.4f}")

    # Test 1: Direct & Maximized Overlap / Mismatch (SXS vs DANSur)
    ifo_h1 = ifos[0]
    h_sxs_fd = signals["H1"]
    mask_f = (ifo_h1.frequency_array >= ifo_h1.minimum_frequency) & (ifo_h1.frequency_array <= ifo_h1.maximum_frequency)
    psd = ifo_h1.power_spectral_density_array[mask_f]
    df = ifo_h1.frequency_array[1] - ifo_h1.frequency_array[0]

    def _inner_prod(a, b):
        return 4.0 * df * np.real(np.sum(a[mask_f] * np.conj(b[mask_f]) / psd))

    norm_sxs = _inner_prod(h_sxs_fd, h_sxs_fd)
    norm_dan = _inner_prod(h_dansur_fd, h_dansur_fd)
    direct_overlap = _inner_prod(h_sxs_fd, h_dansur_fd) / np.sqrt(norm_sxs * norm_dan)

    integrand_val = 4.0 * df * (h_sxs_fd[mask_f] * np.conj(h_dansur_fd[mask_f]) / psd)
    N_full_val = (len(h_sxs_fd) - 1) * 2
    full_pos_val = np.zeros(N_full_val, dtype=np.complex128)
    full_pos_val[:len(mask_f)][mask_f] = integrand_val
    corr_val = np.fft.ifft(full_pos_val) * N_full_val
    idx_max = int(np.argmax(np.abs(corr_val)))
    max_overlap = float(np.abs(corr_val[idx_max]) / np.sqrt(norm_sxs * norm_dan))
    dt_res = (idx_max if idx_max < N_full_val // 2 else idx_max - N_full_val) / fs
    dphi_res = float(np.angle(corr_val[idx_max]))

    print(f"  Direct Overlap          : {direct_overlap:.6f} (Mismatch: {1.0 - direct_overlap:.2e})")
    print(f"  Max Overlap (t_c, phi)  : {max_overlap:.6f} (dt: {dt_res*1000:.3f} ms, dphi: {dphi_res:+.4f} rad)")

    min_overlap = _get_float(dicc, "min_overlap", 0.995)
    if max_overlap < min_overlap:
        raise ValueError(
            f"Simulation {sxs_id} rejected: Max Overlap = {max_overlap:.6f} < {min_overlap:.4f} "
            f"(would cause systematic bias in mass_ratio and spins)."
        )

    # Test 2: Mode-by-mode amplitude & phase consistency at t = -0.02 s
    t_dan = np.asarray(domain_dansur)
    t_sxs_aligned = t_uniform + dt_align
    for (l_v, m_v) in [(2, 2), (3, 3), (2, 1), (4, 4)]:
        if (l_v, m_v) in h_modes_dan and any(tuple(h.LM[idx_m]) == (l_v, m_v) for idx_m in keep_indices):
            idx_k = [k for k, idx_m in enumerate(keep_indices) if tuple(h.LM[idx_m]) == (l_v, m_v)][0]
            parity = float(best_parity_sign ** abs(int(m_v)))
            hlm_sxs_aligned = factor * parity * (modes_u_r[idx_k] + 1j * modes_u_i[idx_k]) * np.exp(-1j * m_v * best_dphi)
            idx_sxs_t = int(np.argmin(np.abs(t_sxs_aligned - (-0.02))))
            idx_dan_t = int(np.argmin(np.abs(t_dan - (-0.02))))
            amp_ratio = np.abs(hlm_sxs_aligned[idx_sxs_t]) / (np.abs(h_modes_dan[(l_v, m_v)][idx_dan_t]) + 1e-35)
            dphi_mode = float(np.angle(hlm_sxs_aligned[idx_sxs_t] * np.conj(h_modes_dan[(l_v, m_v)][idx_dan_t])))
            print(f"  Mode ({l_v},{m_v})              -> amp_ratio: {amp_ratio:.4f} | dphi: {dphi_mode:+.4f} rad")

    # ---------------------------------------------------------
    # 8. Save diagnostic waveform plots
    # ---------------------------------------------------------
    times = start_time + np.arange(N) / fs
    t_rel = times - geocent_time
    plt.figure(figsize=(9, 4))
    plt.plot(t_rel, signal_h1_td, label="H1 strain", lw=1.2)
    plt.xlim(-0.5, 0.06)
    plt.xlabel(r"$t - t_c$ (s)")
    plt.ylabel("Strain")
    plt.title("SXS Signal Projected on H1")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("signal_H1_strain.png", dpi=200)
    plt.close()

    t_sxs_rel = t_seconds - t_seconds[np.argmax(np.abs(h_plus_sxs + 1j * h_cross_sxs))]
    plt.figure(figsize=(9, 4))
    plt.plot(t_sxs_rel, h_plus_sxs, label=r"$h_+$ (SXS)", lw=1.2)
    plt.plot(t_sxs_rel, h_cross_sxs, label=r"$h_\times$ (SXS)", lw=1.2, alpha=0.85)
    plt.xlim(-0.5, 0.06)
    plt.xlabel(r"$t - t_{\mathrm{merger}}$ (s)")
    plt.ylabel("Strain")
    plt.title("Polarizations of the gravitational wave (SXS)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("signal_SXS_modes.png", dpi=200)
    plt.close()

    return ifos