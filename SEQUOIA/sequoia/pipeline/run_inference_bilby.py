import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["NUMBA_THREADING_LAYER"] = "workqueue"
os.environ["NUMBA_NUM_THREADS"] = "1"
try:
    import threadpoolctl
    threadpoolctl.threadpool_limits(limits=1, user_api='blas')
    threadpoolctl.threadpool_limits(limits=1, user_api='openmp')
    # print("-> Hilos limitados con threadpoolctl con éxito.")
except ImportError:
    pass 
import bilby
import os
import numpy as np
import time
import json

from bilby.gw.likelihood import GravitationalWaveTransient
from matplotlib import pyplot as plt
from gwpy.timeseries import TimeSeries

#Load differente waveform generators
from sequoia.waveforms.dansur.waveform_generator import waveform_generator
from sequoia.waveforms.nrybsur.nrh_waveform_generator import nrh_waveform_generator
from sequoia.waveforms.nrsur7dq4.nrsur_waveform_generator import nrsur_waveform_generator


# ---------------------------------------------------------
# 8. Run inference
# ---------------------------------------------------------


def run_inference(ifos, dicc,targ_keys, priors, outdir,folder):
    """
    Run Bayesian parameter estimation for a gravitational-wave event.

    The function selects the waveform model specified in ``dicc``, creates
    the corresponding Bilby waveform generator, constructs a
    ``GravitationalWaveTransient`` likelihood, and runs parameter estimation
    with the Nessai sampler. It then saves the execution time, produces a
    corner plot, reconstructs the maximum-likelihood waveform, saves the
    reconstructed waveform plot, and computes the corresponding waveform
    difference diagnostics.

    INPUTS:
    ------
    
    -ifos : bilby.gw.detector.InterferometerList
        List of interferometers containing the detector data, power spectral
        densities, sampling information, and detector strain data for the
        event.

    -dicc : The created dicctionay in config.py 

    -targ_keys : iterable
        Names of the source parameters used by the waveform generator.
        These parameters define the subset of parameters passed from the
        Bilby parameter dictionary to the waveform model.

    -priors : bilby.core.prior.PriorDict
        Prior distributions for the parameters sampled during the Bayesian
        inference.

    -outdir :
        Directory in which Bilby output files, sampler results, execution
        time information, and diagnostic plots are saved.

    -folder : str
        Identifier of the event or simulation folder. It is used to label
        the reconstructed waveform plot.

    Returns
    -------
    -Result:

        Bilby result object containing the posterior samples, likelihood
        evaluations, evidence estimates, and sampler metadata.

    - A JSON file containing the posterior medians and the lower and upper bounds of their 68% credible intervals, 
    together with the injection parameters of the SXS event.   
    """
    # ---------------------------------------------------------
    # 8.1 Checks if DANSur or an approximant will be used
    # ---------------------------------------------------------

    print('=' * 50)
    print('Loading waveform generator...')
    print(dicc["surrogate_model"].lower())

    if dicc["surrogate_model"].lower() in 'dansur':
        generated_waveform = waveform_generator(dicc,targ_keys)

    elif dicc["surrogate_model"].lower() in 'NRHybSur3dq8'.lower() :
        generated_waveform = nrh_waveform_generator(dicc,targ_keys)


    elif dicc["surrogate_model"].lower() in 'NRSur7dq4'.lower() :        
        generated_waveform = nrsur_waveform_generator(dicc,targ_keys)

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
                label= dicc["label"],
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
    title =     dicc["label"] + "_execution_time.json"
    with open(os.path.join(outdir, title), "w") as f:
            json.dump(timing, f, indent=4)

    print(f"Execution time: {execution_time / 60:.2f} min")

    # ---------------------------------------------------------
    # 8.4 Plots reconstructed waveform and corner plot
    # ---------------------------------------------------------

    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')
    title = "Reconstructed waveform for event " + folder +' with_' +     dicc["label"]

    name = folder +' ' +     dicc["label"]+ '_reconstructed_waveform.png'
    plt.title(title, fontsize=11)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, name ))

    result.plot_corner()
    return result



