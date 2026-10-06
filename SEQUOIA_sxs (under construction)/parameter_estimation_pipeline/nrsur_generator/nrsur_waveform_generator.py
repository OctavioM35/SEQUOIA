import numpy as np
import bilby
from scripts.surrogate.sur_utils import DANSur
from parameter_estimation_pipeline.nrsur_generator.nrsur_convert import nrsur_convert
import matplotlib.pyplot as plt
import time
from scipy.interpolate import CubicSpline

from gwpy.frequencyseries import FrequencySeries
import gwsurrogate as gws
nrsur = gws.LoadSurrogate('NRSur7dq4')

def make_my_gen_func(dicc):


    def my_gen_func(times, **kwargs):
        
        # --------------------------------------------------
        # Conversión de parámetros
        # --------------------------------------------------
        converted_params = (
            bilby.gw.conversion
            .convert_to_lal_binary_black_hole_parameters(kwargs)[0]
        )



        out = nrsur_convert(times, **converted_params)


        
        domain, h, _ = nrsur(
            q=out["q"],
            chiA0=out["chiA0"],
            chiB0=out["chiB0"],
            M=out["M"],
            dist_mpc=out["dist_mpc"],
            f_low=dicc["waveform"]["f_low"],
            f_ref=dicc["waveform"]["f_ref"],
            inclination=out["inclination"],
            dt=1 / dicc["sampling-frequency"],
            phi_ref=out["phi_ref"],
            ellMax = 2,
            units="mks",
        )


        domain = np.asarray(domain)
        h = np.asarray(h).squeeze()

        # --------------------------------------------------
        # Centrar merger en t=0
        # --------------------------------------------------
        domain = domain - domain[-1]
        times_rel = times - times[-1]

        # --------------------------------------------------
        # Interpolación
        # --------------------------------------------------
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

        # Evitar np.nan_to_num sobre arrays innecesariamente
        h_plus = np.where(np.isfinite(h_plus), h_plus, 0.0)
        h_cross = np.where(np.isfinite(h_cross), h_cross, 0.0)

        return {
            "plus": h_plus,
            "cross": h_cross,
        }

    return my_gen_func

def nrsur_waveform_generator(dicc, targ_keys):

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