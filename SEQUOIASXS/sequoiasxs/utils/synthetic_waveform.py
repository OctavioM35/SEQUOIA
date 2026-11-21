import numpy as np
import bilby
import sxs
from scipy.interpolate import CubicSpline
from bilby.gw.detector import InterferometerList
from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
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
    time convention, domain, taper, and Bilby's frequency-domain projection.
    """
    # ---------------------------------------------------------
    # 1. Configuration & defaults from dicc
    # ---------------------------------------------------------
    total_mass = _get_float(dicc, "total_mass", 70.0)                         # M_sun
    luminosity_distance_mpc = _get_float(dicc, "luminosity_distance", 150.0)  # Mpc
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
    except Exception:
        sim_obj = sxs.load(sxs_id)

    h = sim_obj.h if hasattr(sim_obj, "h") else sim_obj

    m1_ref = float(h.metadata.get("reference_mass1"))
    m2_ref = float(h.metadata.get("reference_mass2"))
    if m1_ref is None or m2_ref is None:
        raise ValueError(f"Missing reference masses in SXS metadata for {sxs_id}.")

    spin1_ref = h.metadata.get("reference_dimensionless_spin1", [0.0, 0.0, 0.0])
    spin2_ref = h.metadata.get("reference_dimensionless_spin2", [0.0, 0.0, 0.0])

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
    # 3. Filter modes, match DANSur domain, interpolate & taper uniformly
    # ---------------------------------------------------------
    # 3.1 Query exact DANSur innate time domain (in seconds relative to merger)
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
    domain_dansur, _, _ = sur(**dansur_query)
    t_start_dansur_s = float(domain_dansur[0])
    t_end_dansur_s = float(domain_dansur[-1])

    # 3.2 Filter exact (l, +-m) modes in SXS
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

    # 3.3 Align t=0 with the L2 norm across active modes
    t_all = np.asarray(h.t)
    mode_norm = np.linalg.norm(h_data, axis=1)
    t_peak = t_all[np.argmax(mode_norm)]

    relaxation = float(h.metadata.get("relaxation_time", t_all[0]))
    t_start_sxs_s = (relaxation - t_peak) * time_conversion

    # Start time is the later of DANSur start or SXS relaxation time
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
    # 5. Interpolate onto uniform linear grid, taper start, then wrap cyclically
    # ---------------------------------------------------------
    t_src = t_seconds_all[valid_sxs]
    # Rejilla lineal uniforme cronológica (de inicio de inspiral a fin de ringdown)
    t_uniform = np.arange(t_src[0], t_src[-1], 1.0 / fs)
    n_u = len(t_uniform)
    n_active = len(keep_indices)

    modes_u_r = np.zeros((n_active, n_u), dtype=np.float64)
    modes_u_i = np.zeros((n_active, n_u), dtype=np.float64)

    for k, idx_mode in enumerate(keep_indices):
        hlm_src = h_data[valid_sxs, idx_mode]
        modes_u_r[k] = np.interp(t_uniform, t_src, hlm_src.real)
        modes_u_i[k] = np.interp(t_uniform, t_src, hlm_src.imag)

    # Aplicar el Planck taper en el inicio cronológico (t_uniform[0], inicio del inspiral)
    modes_u_r, modes_u_i = _numba_taper_modes(modes_u_r, modes_u_i)

    # Sumar modos con _sYlm (probamos primero sin offset y con -0.2755)
    from gwsurrogate.surrogate import _sYlm
    phi_eval = float(-phi - 0.2755)
    factor = (G * total_mass * M_sun) / (c**2 * luminosity_distance_mpc * Mpc)

    hp_uniform = np.zeros(n_u, dtype=np.float64)
    hc_uniform = np.zeros(n_u, dtype=np.float64)

    for k, idx_mode in enumerate(keep_indices):
        l_val, m_val = h.LM[idx_mode]
        ylm = _sYlm(-2, int(l_val), int(m_val), float(theta), phi_eval)
        hr = ylm.real * modes_u_r[k] - ylm.imag * modes_u_i[k]
        hi = ylm.real * modes_u_i[k] + ylm.imag * modes_u_r[k]
        hp_uniform += factor * hr
        hc_uniform += -factor * hi

    # Calcular f_min justo después del taper en la malla uniforme
    idx_22_local = [k for k, idx_m in enumerate(keep_indices) if tuple(h.LM[idx_m]) == (2, 2)][0]
    h22_u = modes_u_r[idx_22_local] + 1j * modes_u_i[idx_22_local]
    amp_22 = np.abs(h22_u)
    taper_end = np.where(amp_22 >= 0.95 * np.max(amp_22[:max(10, n_u // 4)]))[0]
    idx_fstart = int(taper_end[0]) + 5 if len(taper_end) > 0 else int(0.1 * n_u)
    idx_fstart = min(idx_fstart, n_u - 2)

    phase_22 = np.unwrap(np.angle(h22_u))
    dphi_dt = np.gradient(phase_22, t_uniform)
    f_start_trimmed = abs(float(dphi_dt[idx_fstart])) / (2.0 * np.pi)
    f_min = float(max(20.0, f_start_trimmed * 1.15))

    dicc.setdefault("waveform", {})
    dicc["waveform"]["minimum-frequency"] = float(f_min)
    dicc["waveform"]["reference-frequency"] = float(f_ref)
    dicc["waveform"]["f_low"] = float(f_min)
    dicc["waveform"]["f_ref"] = float(f_ref)

    # Colocar en el buffer cíclico tau_cyclic (t >= 0 al inicio, t < 0 al final)
    t_elapsed = np.arange(N) / fs
    t_post_max = t_ringdown + 0.01
    tau_cyclic = np.where(t_elapsed <= t_post_max, t_elapsed, t_elapsed - duration)

    hp_cyclic = np.interp(tau_cyclic, t_uniform, hp_uniform, left=0.0, right=0.0)
    hc_cyclic = np.interp(tau_cyclic, t_uniform, hc_uniform, left=0.0, right=0.0)

    # Guardar para diagnósticos
    h_plus_sxs = hp_uniform
    h_cross_sxs = hc_uniform
    t_seconds = t_uniform

    # Transformar al dominio de la frecuencia con nfft
    hp_freq, _ = bilby.core.utils.nfft(hp_cyclic, fs)
    hc_freq, _ = bilby.core.utils.nfft(hc_cyclic, fs)
    polarizations = {"plus": hp_freq, "cross": hc_freq}

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
    # ---------------------------------------------------------
    # 6. Initialize interferometers & inject in frequency domain
    # ---------------------------------------------------------
    freqs = np.linspace(0, fs / 2, N // 2 + 1)
    psd_array = np.ones_like(freqs) * psd_value

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
        ifo.set_strain_data_from_power_spectral_density(
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
    # 7. Store injection dictionary & diagnostics
    # ---------------------------------------------------------
    signal_h1_td = bilby.core.utils.infft(signals["H1"], fs)
    dicc["signal"] = [signal_h1_td, h_cross_sxs, h_plus_sxs]
    dicc["inyection"] = injection_parameters
    dicc["chirp_mass"] = float(chirp_mass)

    snr_dict = {}
    for ifo in ifos:
        snr_sq = ifo.optimal_snr_squared(signals[ifo.name])
        snr_dict[ifo.name] = float(np.sqrt(np.real(snr_sq)))
        print(f"  {ifo.name} optimal SNR: {snr_dict[ifo.name]:.4f}")
    snr_net = np.sqrt(sum(s**2 for s in snr_dict.values()))
    print(f"  Network optimal SNR: {snr_net:.4f}")

    times = start_time + np.arange(N) / fs

    # 1. Gráfica de la señal proyectada en H1
    t_rel = times - geocent_time
    plt.figure(figsize=(9, 4))
    plt.plot(t_rel, signal_h1_td, label="H1 strain", lw=1.2)
    plt.xlim(-0.5, 0.06)
    plt.xlabel(r"$t - t_c$ (s)")
    plt.ylabel("Strain")
    plt.title("Señal SXS proyectada en H1")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("signal_H1_strain.png", dpi=200)
    plt.close()

    # 2. Gráfica de las polarizaciones puras de SXS
    t_sxs_rel = t_seconds - t_seconds[np.argmax(np.abs(h_plus_sxs + 1j * h_cross_sxs))]
    plt.figure(figsize=(9, 4))
    plt.plot(t_sxs_rel, h_plus_sxs, label=r"$h_+$ (SXS)", lw=1.2)
    plt.plot(t_sxs_rel, h_cross_sxs, label=r"$h_\times$ (SXS)", lw=1.2, alpha=0.85)
    plt.xlim(-0.5, 0.06)
    plt.xlabel(r"$t - t_{\mathrm{merger}}$ (s)")
    plt.ylabel("Strain")
    plt.title("Polarizaciones de la onda gravitacional (SXS)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("signal_SXS_modes.png", dpi=200)
    plt.close()
    # # 1. Crear el generador de DANSur y evaluar en los parámetros exactos de inyección
    inj_params = dicc["inyection"].copy()
    wf_gen = waveform_generator(dicc, list(inj_params.keys()))

    # Obtener polarizaciones en frecuencia de DANSur y proyectar en H1
    pols_dansur = wf_gen.frequency_domain_strain(inj_params)
    ifo_h1 = ifos[0]  # H1
    h_dansur_fd = ifo_h1.get_detector_response(pols_dansur, inj_params)

    # 2. Señal SXS exacta en frecuencia (sin ruido)
    # Reconstruimos polarizaciones SXS desde el buffer o guardamos signals["H1"]
    h_sxs_fd = bilby.core.utils.nfft(dicc["signal"][0], dicc["sampling-frequency"])[0]

    # 3. Calcular el Overlap y el Mismatch en la banda [f_min, f_max]
    mask_f = (ifo_h1.frequency_array >= ifo_h1.minimum_frequency) & (ifo_h1.frequency_array <= ifo_h1.maximum_frequency)
    psd = ifo_h1.power_spectral_density_array[mask_f]
    df = ifo_h1.frequency_array[1] - ifo_h1.frequency_array[0]

    def inner_product(a, b):
        return 4.0 * df * np.real(np.sum(a[mask_f] * np.conj(b[mask_f]) / psd))

    overlap = inner_product(h_sxs_fd, h_dansur_fd) / np.sqrt(
        inner_product(h_sxs_fd, h_sxs_fd) * inner_product(h_dansur_fd, h_dansur_fd)
    )
    mismatch = 1.0 - overlap

    print(f"Overlap (SXS vs DANSur):  {overlap:.6f}")
    print(f"Mismatch (1 - Overlap):   {mismatch:.2e}")

    def complex_max_overlap(a_fd, b_fd, mask_f, psd, df, fs):
        # Producto interno complejo por frecuencia
        integrand = 4.0 * df * (a_fd[mask_f] * np.conj(b_fd[mask_f]) / psd)
        N_full = (len(a_fd) - 1) * 2
        full_pos = np.zeros(N_full, dtype=np.complex128)
        full_pos[:len(mask_f)][mask_f] = integrand

        # IFFT compleja: el módulo da el overlap maximizado sobre t_c y phase a la vez
        corr = np.fft.ifft(full_pos) * N_full
        norm = np.sqrt(
            4.0 * df * np.real(np.sum(np.abs(a_fd[mask_f])**2 / psd)) *
            4.0 * df * np.real(np.sum(np.abs(b_fd[mask_f])**2 / psd))
        )
        idx = np.argmax(np.abs(corr))
        dt = (idx if idx < N_full // 2 else idx - N_full) / fs
        dphi = np.angle(corr[idx])
        return float(np.abs(corr[idx]) / norm), dt, dphi


    # --- Llamada al diagnóstico ---
    ifo_h1 = ifos[0]
    mask_f = (ifo_h1.frequency_array >= ifo_h1.minimum_frequency) & (ifo_h1.frequency_array <= ifo_h1.maximum_frequency)
    psd = ifo_h1.power_spectral_density_array[mask_f]
    df = ifo_h1.frequency_array[1] - ifo_h1.frequency_array[0]

    # Señal SXS vs señal DANSur en frecuencia (sin ruido)
    h_sxs_fd = signals["H1"]  # o bilby.core.utils.nfft(dicc["signal"][0], fs)[0]
    wf_gen = waveform_generator(dicc, list(injection_parameters.keys()))
    pols_dansur = wf_gen.frequency_domain_strain(injection_parameters)
    h_dansur_fd = ifo_h1.get_detector_response(pols_dansur, injection_parameters)

    ov_true, dt_true, dphi_true = complex_max_overlap(h_sxs_fd, h_dansur_fd, mask_f, psd, df, fs)
    print(f"Max Overlap REAL (t_c + phase): {ov_true:.6f} | dt: {dt_true*1000:.3f} ms | dphi: {dphi_true:.4f} rad")
    # raise Exception("Detener ejecución para inspección de resultados. Ajustar t_c y phi según corresponda.")
    return ifos