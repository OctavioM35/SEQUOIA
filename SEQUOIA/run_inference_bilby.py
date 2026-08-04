
import bilby
from bilby.gw.likelihood import GravitationalWaveTransient
from waveform_generator import  waveform_generator
from matplotlib import pyplot as plt
import os
from gwpy.timeseries import TimeSeries
from approximant_wv_generator import approximant_generator
import sys
import h5py
import numpy as np
import time


def test_likelihood(likelihood, ifos, waveform_generator):

    print("\n" + "=" * 70)
    print("TEST LIKELIHOOD")
    print("=" * 70)

    # -------------------------------
    # Noise
    # -------------------------------

    noise = likelihood.noise_log_likelihood()

    print("Noise logL:", noise)
    print("Finite:", np.isfinite(noise))

    # -------------------------------
    # Parámetros
    # -------------------------------

    parameters = likelihood.parameters.copy()

    print("\nNúmero de parámetros:")
    print(len(parameters))

    # -------------------------------
    # Waveform
    # -------------------------------

    try:

        waveform = (
            waveform_generator
            .frequency_domain_strain(parameters)
        )

        print("\nWaveform OK")

        for key, value in waveform.items():

            print(
                key,
                "size =", len(value),
                "NaN =", np.isnan(value).any(),
                "Inf =", np.isinf(value).any()
            )

    except Exception as e:

        print("\nERROR WAVEFORM:")
        print(e)

    # -------------------------------
    # Likelihood
    # -------------------------------

    try:

        ratio = likelihood.log_likelihood_ratio()

        total = likelihood.log_likelihood()

        print("\nLog likelihood ratio:", ratio)
        print("Log likelihood:", total)

        print(
            "Consistency:",
            total - (noise + ratio)
        )

    except Exception as e:

        print("\nERROR LIKELIHOOD:")
        print(e)

    # -------------------------------
    # Interferómetros
    # -------------------------------

    print("\nINTERFERÓMETROS")

    for ifo in ifos:

        print(
            ifo.name,
            "duration =", ifo.duration,
            "fs =", ifo.sampling_frequency,
            "fmin =", ifo.minimum_frequency,
            "fmax =", ifo.maximum_frequency
        )







def run_inference(ifos, dicc,targ_keys, priors, npoints, outdir,resume,approximant, zenodo_file):

    if approximant== False:
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"

    elif approximant == True:
        generated_waveform = approximant_generator(dicc,targ_keys)
        label = "IMR"

    else:
        print('Approximant (in config.py) must either be True or False.')
        sys.exit(1)


    likelihood = GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=dicc["phase-marginalization"],
        time_marginalization=dicc["time-marginalization"],
        jitter_time=dicc["jitter-time"]
        # calibration_marginalization= dicc["calibration-marginalization"]
    )

        # ==========================================================
    # MONITOR DE LIKELIHOOD
    # ==========================================================

    original_log_likelihood = likelihood.log_likelihood

    likelihood_counter = {
        "n": 0
    }

    likelihood_start = time.time()


    def monitored_log_likelihood():

        likelihood_counter["n"] += 1

        n = likelihood_counter["n"]

        t0 = time.time()

        # ------------------------------------------------------
        # Evaluar likelihood original
        # ------------------------------------------------------

        log_likelihood_value = original_log_likelihood()

        elapsed = time.time() - t0
        total_elapsed = time.time() - likelihood_start


        # ------------------------------------------------------
        # Mostrar progreso
        # ------------------------------------------------------

        print(
            f"[LIKELIHOOD] "
            f"Evaluation #{n} | "
            f"time = {elapsed:.3f} s | "
            f"total = {total_elapsed:.1f} s | "
            f"logL = {log_likelihood_value}",
            flush=True
        )

        return log_likelihood_value


    likelihood.log_likelihood = monitored_log_likelihood


    # ==========================================================
    # Información inicial
    # ==========================================================

    print("\n==============================================")
    print("Starting Bilby / Nessai")
    print("==============================================")

    print(f"Sampler       : nessai")
    print(f"Waveform      : {label}")
    print(f"npoints       : {npoints}")
    print(f"naccept       : 10")
    print(f"npool         : 1")
    print(f"resume        : {resume}")
    print("==============================================\n")


    result=bilby.run_sampler(
                likelihood=likelihood,
                priors=priors,
                sampler="nessai",
                use_ratio=False,
                flow_class="gwflowproposal",
                npoints=npoints,
                seed=1234,
                naccept=10,
                resume=resume,
                outdir=outdir,
                label=label,
                npool=1,
                stopping=1
            )
    
    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    plt.title("Reconstructed waveform", fontsize=14)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, 'reconstructed_waveform.png' ))
    return result



