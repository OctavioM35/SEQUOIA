import numpy as np
import bilby
import matplotlib.pyplot as plt
import time
import gwsurrogate as gws


from scipy.interpolate import CubicSpline
from gwpy.frequencyseries import FrequencySeries
from sequoiasxs.waveforms.nrsur7dq4.nrsur_convert import nrsur_convert

nrsur = gws.LoadSurrogate('NRSur7dq4')

def make_my_gen_func(dicc):
    """
    Build the time-domain source model used by the Bilby waveform generator.

    The function reads the low and reference frequencies from the
    configuration dictionary and returns a closure (my_gen_func) that
    evaluates the NRHybSur3dq8 surrogate for a given set of sampled
    parameters and resamples the result onto the time array of the data.

    INPUTS:
    -------

    -dicc: Configuration dictionary. The keys dicc["waveform"]["f_low"] and
    dicc["waveform"]["f_ref"] (Hz) are used.

    OUTPUTS:
    --------

    -my_gen_func: Function with signature my_gen_func(times, **kwargs) that
    returns a dictionary {"plus": h_plus, "cross": h_cross} with the two
    polarizations evaluated at "times". Specifically, it:
        1. Converts the sampled parameters to the LAL binary black hole
           convention (bilby.gw.conversion).
        2. Maps them to the surrogate arguments with nnsur_convert.
        3. Evaluates the surrogate with aligned spins [0, 0, chi] and
           the f_low and f_ref of the configuration, in SI units.
        4. Shifts the surrogate time axis so that its last sample coincides
           with the last sample of "times".
        5. Interpolates the real and imaginary parts of the strain with
           cubic splines. Times outside the surrogate domain are set to zero.
    """

    def my_gen_func(times, **kwargs):

        converted_params = (
            bilby.gw.conversion
            .convert_to_lal_binary_black_hole_parameters(kwargs)[0]
        )
        out = nrsur_convert(times, **converted_params)
        chiA0 = np.asarray(out["chiA0"], dtype=float)
        chiB0 = np.asarray(out["chiB0"], dtype=float)

        chiA0 = [0.0, 0.0, chiA0[2]]
        chiB0 = [0.0, 0.0, chiB0[2]]

        domain, h, _ = nrsur(
            q=out["q"],
            chiA0=chiA0,
            chiB0=chiB0,
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

def nrsur_waveform_generator(dicc, targ_keys):
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