import numpy as np
import bilby
from scripts.surrogate.sur_utils import DANSur
from waveform_conversion import nnsur_convert
from gwpy.frequencyseries import FrequencySeries

nnsur = DANSur(device="cpu")


def my_gen_func(times, **kwargs):

    converted = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(kwargs)[0]
    out = nnsur_convert(times, **converted)
    domain, h, _ = nnsur(
        q=out["q"],
        chiA0=[0, 0, out["chiA0"]],
        chiB0=[0, 0, out["chiB0"]],
        M=out["M"],
        dist_mpc=out["dist_mpc"],
        times=times,
        f_low=0,
        inclination=out["inclination"],
        phi_ref=out["phi_ref"],
        units="mks"
    )

    h = np.squeeze(np.asarray(h))

    h_plus = np.real(h)
    h_cross = np.imag(h)

    return {"plus": h_plus, "cross": h_cross}

def waveform_generator(duration, sampling_frequency,targ_keys):
        waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            duration=duration,
            sampling_frequency=sampling_frequency,
            # time_domain_source_model=mymodel,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, 
            time_domain_source_model=my_gen_func,
            frequency_domain_source_model=None
            # start_time=waveform_generator2.start_time

        )        
        waveform_generator.source_parameter_keys = targ_keys
        return waveform_generator