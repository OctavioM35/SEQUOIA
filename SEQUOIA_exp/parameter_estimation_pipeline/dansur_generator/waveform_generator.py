import numpy as np
import bilby
from scripts.surrogate.sur_utils import DANSur
from parameter_estimation_pipeline.dansur_generator.waveform_conversion import nnsur_convert
from gwpy.frequencyseries import FrequencySeries
import gwsurrogate as gws
import matplotlib.pyplot as plt

nnsur = DANSur(device="cpu")
import time


def make_my_gen_func(dicc):
    def my_gen_func(times, **kwargs):
 
                converted_params = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(kwargs)[0]


                out = nnsur_convert(times, **converted_params)
                # t0 = time.perf_counter()

                # raise Exception
                domain, h, _ = nnsur(
                    q=out['q'],
                    chiA0=[0, 0, out['chiA0']],
                    chiB0=[0, 0, out['chiB0']],
                    M=out['M'],
                    dist_mpc=out['dist_mpc'],
                    times=times,
                    f_low=dicc["waveform"]["f_low"],
                    inclination=out['inclination'],
                    f_ref=dicc["waveform"]["f_ref"],
                    phi_ref=out['phi_ref'],
                    units='mks'
                )
                # t_eval = time.perf_counter() - t0
                # print(f"Tiempo de evaluación dansur: {t_eval:.6f} s")
                h = np.squeeze(np.asarray(h))
                h_plus = np.real(h)
                h_cross = np.imag(h)
                # print(h)
                # print(len(h))
                # print(type(h))
                # print(len(h))

                # plt.figure(figsize=(10, 5))
                # plt.plot(domain, h_plus, label="h_plus")
                # plt.plot(domain, h_cross, label="h_cross")
                # plt.xlabel("Time [s]")
                # plt.ylabel("Strain")
                # plt.legend()
                # plt.grid()
                # plt.tight_layout()
                # plt.show()
                # plt.savefig('prueba dansur')
                return {"plus": h_plus, "cross": h_cross}
    return my_gen_func

def waveform_generator(dicc,targ_keys):

        waveform_generator = bilby.gw.waveform_generator.WaveformGenerator(
            duration=float(dicc["duration"]),
            sampling_frequency=float(dicc["sampling-frequency"]),
            # time_domain_source_model=mymodel,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, 
            time_domain_source_model=make_my_gen_func(dicc),
            frequency_domain_source_model=None,
            # waveform_arguments=waveform_arguments,
            start_time = dicc["start_time"]
        )        
        waveform_generator.source_parameter_keys = set(targ_keys)

        return waveform_generator