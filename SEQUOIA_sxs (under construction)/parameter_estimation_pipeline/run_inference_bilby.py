import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"


try:
    import threadpoolctl
    threadpoolctl.threadpool_limits(limits=1, user_api='blas')
    threadpoolctl.threadpool_limits(limits=1, user_api='openmp')
    print("-> Hilos limitados con threadpoolctl con éxito.")
except ImportError:
    pass 
import bilby
from bilby.gw.likelihood import GravitationalWaveTransient
from parameter_estimation_pipeline.dansur_generator.waveform_generator import waveform_generator
from matplotlib import pyplot as plt
import os
from gwpy.timeseries import TimeSeries
from parameter_estimation_pipeline.approximant_generator.approximant_wv_generator import approximant_generator
from parameter_estimation_pipeline.nrh_generator.nrh_waveform_generator import nrh_waveform_generator

from parameter_estimation_pipeline.nrsur_generator.nrsur_waveform_generator import nrsur_waveform_generator
import sys
import h5py
import numpy as np
import time
import json
import numpy as np
# ---------------------------------------------------------
# 8. Run inference
# ---------------------------------------------------------




def run_inference(ifos, dicc,targ_keys, priors, outdir,folder):

    # ---------------------------------------------------------
    # 8.1 Checks if DANSur or an approximant will be used
    # ---------------------------------------------------------

    print('=' * 50)
    print('Loading waveform generator...')
    print(dicc["surrogate_model"].lower())

    if dicc["surrogate_model"].lower() in 'dansur':
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"

    elif dicc["surrogate_model"].lower() in 'NRHybSur3dq8'.lower() :
        generated_waveform = nrh_waveform_generator(dicc,targ_keys)
        label = "NRHybSur3dq8"


    elif dicc["surrogate_model"].lower() in 'imr' :
        generated_waveform = approximant_generator(dicc,targ_keys)
        label = "IMR"

    elif dicc["surrogate_model"].lower() in 'NRSur7dq4'.lower() :        
        generated_waveform = nrsur_waveform_generator(dicc,targ_keys)
        label = "NRSur"

    dicc["label"] =label
    print('Waveform generator loaded')
    print('=' * 50)

    # ---------------------------------------------------------
    # 8.2 Runs PE
    # ---------------------------------------------------------
    
    # Force initialization of frequency-domain data and frequency mask

    for ifo in ifos:
            ifo.frequency_domain_strain.shape
            ifo.frequency_mask.shape
            ifo.frequency_array.shape
    class NRSURLikelihood(GravitationalWaveTransient):

        def log_likelihood_ratio(self, parameters):

            try:
                return super().log_likelihood_ratio(parameters)

            except Exception as e:

                if "omega_ref" in str(e):
                    return -np.inf

                raise
    t0 = time.perf_counter()
        
    likelihood = NRSURLikelihood(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=True,
        time_marginalization=True,
        jitter_time=False,

        # calibration_marginalization= dicc["marginalization"]["calibration"]
    )

    dt = time.perf_counter() - t0
    # print('====' * 30)
    # print(f"PID={os.getpid()} | likelihood={dt:.4f} s")


    start_time = time.time()
        # phase_marginalization=dicc["marginalization"]["phase"],
        # time_marginalization=dicc["marginalization"]["time"],
        # jitter_time=dicc["marginalization"]["jitter"],

        # distance_marginalization=False,
        # phase_marginalization=True,
        # time_marginalization=True,
        # jitter_time=False,

    result=bilby.run_sampler(
                likelihood=likelihood,
                priors=priors,
                sampler="nessai",
                use_ratio=False,
                flow_class="gwflowproposal",
                npoints=dicc["inference"]["npoints"],
                resume=dicc["inference"]["resume"],
                # injection_parameters=dicc["injection_parameters"]
                clean=not dicc["inference"]["resume"],
                outdir=outdir,
                label=label,
                npool=14,
                nparallel=14,
                stopping=dicc["inference"]["stopping"],
                sampling_seed=dicc["inference"]["sampling-seed"]
            )

    execution_time = time.time() - start_time

    timing = {
        "execution_time_seconds": execution_time,
        "execution_time_minutes": execution_time / 60,
        "execution_time_hours":  execution_time / 3600,
    }


    # ---------------------------------------------------------
    # 8.3 Saves execution time
    # ---------------------------------------------------------
    title = label + "_execution_time.json"
    with open(os.path.join(outdir, title), "w") as f:
            json.dump(timing, f, indent=4)

    print(f"Execution time: {execution_time / 60:.2f} min")

    # ---------------------------------------------------------
    # 8.5 Plots reconstructed waveform and corner plot
    # ---------------------------------------------------------
    result.plot_corner()

    plt.figure(figsize=(10, 5))
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    generated_wave= generated_waveform.time_domain_strain(maxll_params)['plus']
    dicc["generated_waveform"] = generated_wave
    if len(generated_wave) > len(generated_waveform.time_array):
        plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'][:-1], '--')
    else:
        plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    title = "Reconstructed waveform for event " + folder +' with ' + label
    name = folder + '_reconstructed_waveform.png'
    plt.title(title, fontsize=11)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, name ))
    from scipy.signal import find_peaks


    peaks, properties = find_peaks(generated_wave)

    x_peaks = generated_wave[peaks]
    y = generated_waveform.time_array[peaks]

    diferencias = np.diff(y)
    media_diferencias = np.mean(diferencias)

    print("Diferencias:", diferencias)
    print("Media:", media_diferencias)
    print('maximos', y)    
    return result



