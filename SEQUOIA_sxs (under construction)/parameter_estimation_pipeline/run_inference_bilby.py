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
from parameter_estimation_pipeline.dif import diference

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
    # plt.plot(ifos[0].time_domain_strain)
    # plt.savefig('checkpoint1.png')
    # plt.close()
    # plt.loglog(ifos[0].frequency_array, ifos[0].power_spectral_density_array)
    # plt.savefig('checkpoint2.png')
    # raise Exception
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
                # injection_parameters=dicc["inyection"],
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
    # 8.3 Save execution time
    # ---------------------------------------------------------


    title = label + "_execution_time.json"
    time_outdir = os.path.join(outdir, title)
    with open(time_outdir, "w") as f:
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
    plt.title(title, fontsize=11)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, title ))
    diference(outdir,label,dicc)


    return result



