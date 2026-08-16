#run_inference_bilby
import bilby
from bilby.gw.likelihood import GravitationalWaveTransient
from parameter_estimation_pipeline.waveform_generator import waveform_generator
from matplotlib import pyplot as plt
import os
from gwpy.timeseries import TimeSeries
from parameter_estimation_pipeline.approximant_wv_generator import approximant_generator
import sys
import h5py
import numpy as np
import time
import json

def run_inference(ifos, dicc,targ_keys, priors, outdir, zenodo_file,folder):

    # ---------------------------------------------------------
    # 8.1 Checks if DANSur or an approximant will be used
    # ---------------------------------------------------------

    if dicc["inference"]["approximant"] == False:
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"

    elif dicc["inference"]["approximant"] == True:
        generated_waveform = approximant_generator(dicc,targ_keys)
        label = "IMR"

    else:
        raise ValueError('approximant must either be True or False.')


    # ---------------------------------------------------------
    # 8.2 Inyects synthetic signal (if used)
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

    likelihood = GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=True,
        time_marginalization=True,
        jitter_time=False
        # calibration_marginalization= dicc["calibration-marginalization"]
    )
    start_time = time.time()

    result=bilby.run_sampler(
                likelihood=likelihood,
                priors=priors,
                sampler="nessai",
                use_ratio=False,
                flow_class="gwflowproposal",
                npoints=dicc["inference"]["npoints"],
                resume=dicc["inference"]["resume"],
                clean=not dicc["inference"]["resume"],
                outdir=outdir,
                label=label,
                npool=1,
                # nparallel=2,
                stopping=dicc["inference"]["stopping"],
                sampling_seed=dicc["inference"]["sampling-seed"]
            )

    execution_time = time.time() - start_time

    timing = {
        "execution_time_seconds": execution_time,
        "execution_time_minutes": execution_time / 60,
    }

    with open(os.path.join(outdir, "execution_time.json"), "w") as f:
        json.dump(timing, f, indent=4)

    print(f"Execution time: {execution_time / 60:.2f} min")
    # ---------------------------------------------------------
    # 8.4 Plots reconstructed waveform and corner plot
    # ---------------------------------------------------------

    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    plt.title("Reconstructed waveform", fontsize=14)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, 'reconstructed_waveform.png' ))
    return result



