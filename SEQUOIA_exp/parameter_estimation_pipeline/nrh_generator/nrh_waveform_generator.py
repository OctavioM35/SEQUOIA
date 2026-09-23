import numpy as np
import bilby
from scripts.surrogate.sur_utils import DANSur
from parameter_estimation_pipeline.nrh_generator.nrh_waveform_conversion import nnsur_convert
from parameter_estimation_pipeline.nrh_generator.nrh_wrapper  import NRHybSur3dq8
import matplotlib.pyplot as plt
import time
from scipy.interpolate import CubicSpline


from gwpy.frequencyseries import FrequencySeries
import gwsurrogate as gws
nrh = gws.LoadSurrogate("NRHybSur3dq8")


def make_my_gen_func(dicc):

    f_low = float(dicc["waveform"]["f_low"])
    f_ref = float(dicc["waveform"]["f_ref"])

    def my_gen_func(times, **kwargs):

        converted_params = (
            bilby.gw.conversion
            .convert_to_lal_binary_black_hole_parameters(kwargs)[0]
        )

        out = nnsur_convert(times, **converted_params)

        domain, h, _ = nrh(
            q=out["q"],
            chiA0=[0.0, 0.0, out["chiA0"]],
            chiB0=[0.0, 0.0, out["chiB0"]],
            M=out["M"],
            dist_mpc=out["dist_mpc"],
            f_low=f_low,
            inclination=out["inclination"],
            f_ref=f_ref,
            phi_ref=out["phi_ref"],
            mode_list=[(2, 2)],
            units="mks",
        )

        domain = np.asarray(domain)
        h = np.asarray(h).squeeze()

        domain = domain - domain[-1]
        times_rel = times - times[-1]


        cs_plus = CubicSpline(
            domain,
            h.real,
            extrapolate=False
        )

        cs_cross = CubicSpline(
            domain,
            h.imag,
            extrapolate=False
        )

        h_plus = cs_plus(times_rel)
        h_cross = cs_cross(times_rel)

        h_plus = np.where(np.isfinite(h_plus), h_plus, 0.0)
        h_cross = np.where(np.isfinite(h_cross), h_cross, 0.0)

        return {
            "plus": h_plus,
            "cross": h_cross,
        }

    return my_gen_func

def nrh_waveform_generator(dicc, targ_keys):

    waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
        duration=float(dicc["duration"]),
        sampling_frequency=float(dicc["sampling-frequency"]),
        parameter_conversion=(
            bilby.gw.conversion
            .convert_to_lal_binary_black_hole_parameters
        ),
        time_domain_source_model=make_my_gen_func(dicc),
        frequency_domain_source_model=None,
        start_time=float(dicc["start_time"]),
    )

    waveform_generator.source_parameter_keys = set(targ_keys)

    return waveform_generator
