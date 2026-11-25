import numba
import numpy as np
import bilby
import torch
torch.set_num_threads(1)
torch.set_num_interop_threads(1)
from numpy.typing import NDArray
from gwsurrogate.surrogate import _sYlm
from functools import partial
from typing import Any
from dansur.surrogate.sur_utils import DANSur
from sequoiasxs.waveforms.dansur.waveform_conversion import nnsur_convert
import matplotlib.pyplot as plt

MODEL_PATH = "models/DANSur_HM_jitted.pt"  # path relative to repo root

modes = modes = [(2, 2), (3, 3), (2, 1), (4, 4), (3, 2), (5, 5), (4, 3)]
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
def _numba_taper_modes(modes_real, modes_imag):
    """Taper the start of each time-aligned mode (parallel over modes)."""
    n_modes, n = modes_real.shape
    out_real = np.empty_like(modes_real)
    out_imag = np.empty_like(modes_imag)

    for k in numba.prange(n_modes):
        dr = modes_real[k].copy()
        di = modes_imag[k].copy()
        dr, di = _numba_taper_complex(dr, di)

        for j in range(n):
            out_real[k, j] = dr[j]
            out_imag[k, j] = di[j]

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

    # 1. Evaluar DANSur directamente con inclination y phi_ref en su dominio innato
    out_eval: dict[str, Any] = out.copy()
    out_eval['f_low'] = 0.0
    out_eval['q'] = out_eval['q'] if out_eval['q'] <= 1 else 1 / out_eval['q']
    out_eval['times'] = None

    # DANSur devuelve h = h_+ - i h_x como un array complejo 1D cuando inclination != None
    domain_full, h_innate, _ = model(**out_eval)

    # 2. Aplicar Planck taper en el dominio cronológico innato ANTES de envolver cíclicamente
    h_r_in = np.ascontiguousarray(h_innate.real[None, :], dtype=np.float64)
    h_i_in = np.ascontiguousarray(h_innate.imag[None, :], dtype=np.float64)
    tapered_r, tapered_i = _numba_taper_modes(h_r_in, h_i_in)

    # 3. Envolver cíclicamente e interpolar a times_arr
    t_post_max = float(domain_full[-1]) + 0.01
    t_elapsed = times_arr - times_arr[0]
    times_relative = np.where(t_elapsed <= t_post_max, t_elapsed, t_elapsed - duration_local)

    h_real = np.interp(times_relative, domain_full, tapered_r[0], left=0.0, right=0.0)
    h_imag = np.interp(times_relative, domain_full, tapered_i[0], left=0.0, right=0.0)

    # En convención NR/gwsurrogate: h = h_+  -i h_x
    return {"plus": -h_real, "cross": h_imag}

def make_my_gen_func(dicc):
    """Returns a plain function f(times, **params) as bilby requires."""
    def gen_func(times, **kwargs):
        return model_gen_func(sur, times, **kwargs)
    return gen_func


def waveform_generator(dicc, targ_keys):
    """
        Create the Bilby WaveformGenerator that uses the DANSurHM surrogate.
    
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