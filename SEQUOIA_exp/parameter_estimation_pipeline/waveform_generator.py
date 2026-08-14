import numpy as np
import bilby
from scripts.surrogate.sur_utils import DANSur
from parameter_estimation_pipeline.waveform_conversion import nnsur_convert
from gwpy.frequencyseries import FrequencySeries
import os
# os.environ["BILBY_INCORRECT_PSD_NORMALIZATION"] = "TRUE"
nnsur = DANSur(device="cpu")


def my_gen_func(times, **kwargs):
            converted_params = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(kwargs)[0]
            out = nnsur_convert(times, **converted_params)

            #times = waveform_generator.time_array   REVISAR MAÑANA ESTA LINEA

            domain, h, _ = nnsur(
                q=out['q'],
                chiA0=[0, 0, out['chiA0']],
                chiB0=[0, 0, out['chiB0']],
                M=out['M'],
                dist_mpc=out['dist_mpc'],
                times=times,
                f_low=0,
                inclination=out['inclination'],
                phi_ref=out['phi_ref'],
                units='mks'
            )

            h = np.squeeze(np.asarray(h))

            h_plus = np.real(h)
            h_cross = np.imag(h)

            return {"plus": h_plus, "cross": h_cross}


def waveform_generator(dicc,targ_keys):
        waveform_arguments = dict(
            # minimum_frequency=int(dicc["waveform"]["minimum-frequency"]),
            pn_spin_order=dicc["waveform"]["pn-spin-order"],
            pn_phase_order=dicc["waveform"]["pn-phase-order"],
            pn_tidal_order=dicc["waveform"]["pn-tidal-order"],
            pn_amplitude_order=dicc["waveform"]["pn-amplitude-order"],
            mode_array=dicc["waveform"]["mode-array"],
            catch_waveform_errors=dicc["waveform"]["catch-waveform-errors"],
            reference_frequency=dicc["waveform"]["reference-frequency"],
        )
        waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            duration=float(dicc["duration"]),
            sampling_frequency=float(dicc["sampling-frequency"]),
            # time_domain_source_model=mymodel,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, 
            time_domain_source_model=my_gen_func,
            frequency_domain_source_model=None,
            waveform_arguments=waveform_arguments,
            start_time=dicc["gps_time"]
        )        

        waveform_generator.source_parameter_keys = set(targ_keys)

        return waveform_generator