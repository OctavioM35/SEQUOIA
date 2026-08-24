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


# ---------------------------------------------------------
# 8. Run inference
# ---------------------------------------------------------




def run_inference(ifos, dicc,targ_keys, priors, outdir,folder):

    # ---------------------------------------------------------
    # 8.1 Checks if DANSur or an approximant will be used
    # ---------------------------------------------------------

    print('=' * 50)
    print('Loading waveform generator...')
    if dicc["inference"]["approximant"] == False:
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"
    elif dicc["inference"]["approximant"] == True:
        generated_waveform = approximant_generator(dicc,targ_keys)
        label = "IMR"

    else:
        raise ValueError('approximant must either be True or False.')
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
 
    likelihood = GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=True,
        time_marginalization=True,
        jitter_time=False,

        # calibration_marginalization= dicc["marginalization"]["calibration"]
    )


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
                npool=1,
                # nparallel=2,
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
    if dicc["inference"]["resume"] != True and  dicc["inference"]["skip"] != False :
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
    plt.title(title, fontsize=12)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, name ))

    result.plot_corner()
    return result



