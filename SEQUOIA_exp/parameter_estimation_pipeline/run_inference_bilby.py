import os

# Forzar a OpenBLAS, MKL y la librería científica de C a usar un solo hilo
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

# 2. Si usas librerías que gestionan hilos por debajo, añadimos este seguro técnico:
try:
    import threadpoolctl
    threadpoolctl.threadpool_limits(limits=1, user_api='blas')
    threadpoolctl.threadpool_limits(limits=1, user_api='openmp')
    print("-> Hilos limitados con threadpoolctl con éxito.")
except ImportError:
    pass # Si no tienes la librería instalada, no pasa nada
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
    # 8.2 Inyects synthetic signal (if generate_synthetic_signal == True, see config.py )
    # ---------------------------------------------------------

    if dicc["generate_synthetic_signal"] == True:

            print()
            print("=" * 70)
            print("INJECTING SIGNAL")
            print("=" * 70)

            ifos.inject_signal(
                waveform_generator=generated_waveform,
                parameters=dicc["injection_parameters"]
            )

            print("Injection completed.")

            polarizations = (
                generated_waveform.frequency_domain_strain(
                    dicc["injection_parameters"]
                )
            )

            for ifo in ifos:

                detector_response = (
                    ifo.get_detector_response(
                        polarizations,
                        dicc["injection_parameters"]
                    )
                )


    # ---------------------------------------------------------
    # 8.3 Runs PE
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
                    # print("NRSUR: punto fuera del rango de omega_ref -> logL = -inf")
                    return -np.inf

                # No ocultar otros errores
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
    print('====' * 30)
    print(f"PID={os.getpid()} | likelihood={dt:.4f} s")

    # _original_log_likelihood = likelihood.log_likelihood

    # def debug_log_likelihood(parameters):
    #     t0 = time.time()

    #     print("\n>>> Evaluando likelihood...", flush=True)

    #     try:
    #         logl = _original_log_likelihood(parameters)

    #         elapsed = time.time() - t0

    #         print(
    #             f">>> logL = {logl} | "
    #             f"tiempo = {elapsed:.3f} s | "
    #             f"finite = {np.isfinite(logl)}",
    #             flush=True
    #         )

    #         return logl

    #     except Exception as e:
    #         print(
    #             f">>> ERROR en likelihood después de {time.time()-t0:.3f} s: "
    #             f"{type(e).__name__}: {e}",
    #             flush=True
    #         )
    #         raise

    # likelihood.log_likelihood = debug_log_likelihood

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
    # 8.4 Saves execution time
    # ---------------------------------------------------------
    title = label + "_execution_time.json"
    with open(os.path.join(outdir, title), "w") as f:
            json.dump(timing, f, indent=4)

    print(f"Execution time: {execution_time / 60:.2f} min")

    # ---------------------------------------------------------
    # 8.5 Plots reconstructed waveform and corner plot
    # ---------------------------------------------------------

    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')
    title = "Reconstructed waveform for event " + folder +' with ' + label
    name = folder + '_reconstructed_waveform.png'
    plt.title(title, fontsize=11)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, name ))

    result.plot_corner()
    return result



