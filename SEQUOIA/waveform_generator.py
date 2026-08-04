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
        f_low=11,
        inclination=out["inclination"],
        phi_ref=20,
        units="mks"
    )


    h = np.squeeze(np.asarray(h))

    h_plus = np.real(h)
    h_cross = np.imag(h)
    converted = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(kwargs)[0]

    return {"plus": h_plus, "cross": h_cross}

def waveform_generator(dicc,targ_keys):
        waveform_arguments = dict(
            reference_frequency=dicc["reference-frequency"],
            # minimum_frequency=int(dicc["minimum-frequency"][" waveform"]),
            pn_spin_order=dicc["pn-spin-order"],
            pn_phase_order=dicc["pn-phase-order"],
            pn_tidal_order=dicc["pn-tidal-order"],
            pn_amplitude_order=dicc["pn-amplitude-order"],
            # mode_array=dicc["mode-array"],
            catch_waveform_errors=dicc["catch-waveform-errors"]
        )
        waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            duration=dicc["duration"],
            sampling_frequency=dicc["sampling-frequency"],
            # time_domain_source_model=mymodel,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, 
            time_domain_source_model=my_gen_func,
            frequency_domain_source_model=None,
            waveform_arguments=waveform_arguments,
            start_time=dicc["trigger-time"] - dicc["duration"]
        )        
        waveform_generator.source_parameter_keys = targ_keys
        return waveform_generator