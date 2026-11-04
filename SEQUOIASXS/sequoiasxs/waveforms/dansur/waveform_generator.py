import numba
import numpy as np
import bilby

from numpy.typing import NDArray
from gwsurrogate.surrogate import _sYlm
from functools import partial
from typing import Any
from dansur.surrogate.sur_utils import DANSur
from sequoiasxs.waveforms.dansur.waveform_conversion import nnsur_convert

MODEL_PATH = "models/DANSur_HM_jitted.pt"  # path relative to repo root

modes = [(2, 2), (3, 3), (2, 1), (4, 4), (3, 2), (5, 5), (4, 3)]
sur = DANSur(path=MODEL_PATH, modes_list=modes)


@numba.njit
def _numba_taper_complex(data_real, data_imag):
    """Planck taper at the start of a complex mode (real and imag in one pass).
    Returns (tapered_real, tapered_imag)."""
    length = len(data_real)
    start = -1
    for i in range(length):
        if data_real[i] != 0.0 or data_imag[i] != 0.0:
            start = i
            break
    if start == -1:
        return data_real, data_imag
    end = -1
    for i in range(length - 1, -1, -1):
        if data_real[i] != 0.0 or data_imag[i] != 0.0:
            end = i
            break
    if end - start <= 1:
        return data_real, data_imag
    mid = (start + end) // 2
    flag = 0
    i = start + 1
    n = 0
    LALSIMULATION_RINGING_EXTENT = 19
    amp_i = (data_real[i] ** 2 + data_imag[i] ** 2) ** 0.5
    while flag < 2 and i != mid:
        amp_im1 = (data_real[i - 1] ** 2 + data_imag[i - 1] ** 2) ** 0.5
        amp_ip1 = (data_real[i + 1] ** 2 + data_imag[i + 1] ** 2) ** 0.5
        if amp_i >= amp_im1 and amp_i >= amp_ip1:
            if amp_i == amp_ip1:
                i += 1
                amp_i = amp_ip1
            if i - start > LALSIMULATION_RINGING_EXTENT:
                flag += 1
            n = i - start
        i += 1
        if i < length:
            amp_i = (data_real[i] ** 2 + data_imag[i] ** 2) ** 0.5
    if flag < 2:
        n = mid - start
    if n > 1:
        data_real[start] = 0.0
        data_imag[start] = 0.0
        realN = float(n)
        for j in range(start + 1, start + n - 1):
            realI = float(j - start)
            z = (realN - 1.0) / realI + (realN - 1.0) / (realI - (realN - 1.0))
            sigma = 1.0 / (np.exp(z) + 1.0)
            data_real[j] = data_real[j] * sigma
            data_imag[j] = data_imag[j] * sigma
    return data_real, data_imag


@numba.njit(parallel=True)
def _numba_taper_roll_align_modes(modes_real, modes_imag, roll_by, m_vals, phi22_peak):
    """Fused kernel: taper + cyclic roll left by roll_by + rotation by
    exp(-i * (m/2) * phi22_peak) for every mode (parallel over modes)."""
    n_modes, n = modes_real.shape
    out_real = np.empty_like(modes_real)
    out_imag = np.empty_like(modes_imag)

    for k in numba.prange(n_modes):
        # 1. Taper
        dr = modes_real[k].copy()
        di = modes_imag[k].copy()
        dr, di = _numba_taper_complex(dr, di)

        # 2. Roll left (np.roll(arr, -roll_by))
        roll_by_mod = roll_by % n
        tmp_r = np.empty(n)
        tmp_i = np.empty(n)
        for j in range(n):
            src = (j + roll_by_mod) % n
            tmp_r[j] = dr[src]
            tmp_i[j] = di[src]

        # 3. Phase rotation
        angle = -(m_vals[k] / 2.0) * phi22_peak
        cos_a = np.cos(angle)
        sin_a = np.sin(angle)
        for j in range(n):
            r = tmp_r[j]
            im = tmp_i[j]
            out_real[k, j] = r * cos_a - im * sin_a
            out_imag[k, j] = r * sin_a + im * cos_a

    return out_real, out_imag


@numba.njit(parallel=True)
def _numba_mode_sum_from_precomputed(modes_real, modes_imag, ylm_real, ylm_imag):
    """h = sum_{l,m} Y_{l,m}(theta, phi) * h_{l,m}, parallel over samples.
    The (-m) modes are already in the array, so no symmetry correction.
    Returns (h_real, h_imag) -> (plus, cross)."""
    n_modes, n = modes_real.shape
    h_real = np.zeros(n)
    h_imag = np.zeros(n)

    for j in numba.prange(n):
        hr = 0.0
        hi = 0.0
        for k in range(n_modes):
            ar = ylm_real[k]
            ai = ylm_imag[k]
            mr = modes_real[k, j]
            mi = modes_imag[k, j]
            hr += ar * mr - ai * mi
            hi += ar * mi + ai * mr
        h_real[j] = hr
        h_imag[j] = hi

    return h_real, h_imag


def model_gen_func(model: DANSur, times: NDArray | None = None, **kwargs: Any) -> dict[str, np.ndarray]:
    """Time-domain source model for bilby. Meant to be used via
    ``partial(model_gen_func, model)`` so bilby calls it as ``f(times, **params)``.

    Returns dict with 'plus' and 'cross' (float64 arrays).
    """
    if len(kwargs) == 0:
        raise ValueError('model_gen_func received no parameters')

    if times is None:
        duration_local = 3.0
        times_arr = np.arange(0, duration_local, 1.0 / 1024.0)
    else:
        times_arr = times
        dt = times_arr[1] - times_arr[0]
        duration_local = len(times_arr) * dt

    if 'parameters' in kwargs:
        kwargs = kwargs['parameters']

    converted_params = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(kwargs)[0]
    out = nnsur_convert(times_arr, **converted_params)

    # 1. Generate modes (inclination=None)
    out_no_inc: dict[str, Any] = out.copy()
    out_no_inc['inclination'] = None
    out_no_inc['phi_ref'] = 0.0
    out_no_inc['f_low'] = 0.0

    # convention: 0 < q <= 1
    out_no_inc['q'] = out_no_inc['q'] if out_no_inc['q'] <= 1 else 1 / out_no_inc['q']

    # Evaluate surrogate on its innate domain
    out_no_inc['times'] = None
    domain_full, h_modes_full, _ = model(**out_no_inc)

    # Times relative to a fixed merger epoch (ignores sampled geocent_time
    # to avoid GPS time-scale errors).
    post_trig = domain_full[-1]
    t_merger = float(times_arr[0]) + duration_local - post_trig
    times_relative = times_arr - t_merger

    # Interpolate each mode onto the requested grid
    h_modes: dict[tuple[int, int], np.ndarray] = {}
    for mode, hlm in h_modes_full.items():
        h_modes[mode] = np.interp(times_relative, domain_full, hlm, left=0, right=0)

    # PN relative-phase offsets are already applied by DANSur.__call__.

    # Stack modes into contiguous arrays for Numba
    mode_list_ordered = list(h_modes.keys())
    n_modes = len(mode_list_ordered)
    n_samples = len(times_arr)
    modes_real = np.empty((n_modes, n_samples), dtype=np.float64)
    modes_imag = np.empty((n_modes, n_samples), dtype=np.float64)
    m_vals = np.empty(n_modes, dtype=np.int64)
    for k, (l, m) in enumerate(mode_list_ordered):
        modes_real[k] = h_modes[(l, m)].real
        modes_imag[k] = h_modes[(l, m)].imag
        m_vals[k] = m

    # Index of the merger (closest to t = 0)
    peak_idx = int(np.argmin(np.abs(times_relative)))

    # Taper + cyclic roll, no peak-phase rotation (phi22_peak = 0.0): the
    # DANSur frame is already well defined, and rotating by arg h22(peak)
    # flips the sign of odd-m modes whenever the phase crosses +-pi.
    out_r, out_i = _numba_taper_roll_align_modes(
        modes_real.copy(), modes_imag.copy(), peak_idx, m_vals, 0.0
    )

    # Mode sum with spin-weighted spherical harmonics (LAL convention)
    theta_jn = float(out['inclination'])
    phase = float(out['phi_ref'])
    phi_eval = float(np.pi / 2.0 - phase)
    ylm_r = np.empty(n_modes)
    ylm_i = np.empty(n_modes)
    for k, (l, m) in enumerate(mode_list_ordered):
        y = _sYlm(-2, l, m, theta_jn, phi_eval)
        ylm_r[k] = y.real
        ylm_i[k] = y.imag

    h_real, h_imag = _numba_mode_sum_from_precomputed(out_r, out_i, ylm_r, ylm_i)

    return {"plus": h_real, "cross": h_imag}


def make_my_gen_func(dicc):
    """Returns a plain function f(times, **params) as bilby requires."""
    def gen_func(times, **kwargs):
        return model_gen_func(sur, times, **kwargs)
    return gen_func


def waveform_generator(dicc, targ_keys):
    """
        Create the Bilby WaveformGenerator that uses the NRSur7dq4 surrogate.
    
        The generator works in the time domain, with the duration, sampling
        frequency and start time taken from the configuration dictionary, and
        uses the function built by make_my_gen_func as source model. The
        parameter conversion is the standard LAL binary black hole conversion
        of Bilby.
    
        INPUTS:
        -------
    
        -dicc: Configuration dictionary. The keys "duration", "sampling-frequency"
        and "start_time", plus the entries required by make_my_gen_func, are used.
    
        -targ_keys: Names of the parameters sampled in the inference. They are
        stored as the source parameter keys of the generator.
    
        OUTPUTS:
        --------
    
        -bilby.gw.waveform_generator.WaveformGenerator ready to be passed to the
        Bilby likelihood.
    """
    wf_gen = bilby.gw.waveform_generator.WaveformGenerator(
        duration=float(dicc["duration"]),
        sampling_frequency=float(dicc["sampling-frequency"]),
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
        time_domain_source_model=make_my_gen_func(dicc),
        frequency_domain_source_model=None,
        start_time=dicc["start_time"],
    )
    wf_gen.source_parameter_keys = set(targ_keys)
    return wf_gen