
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

def run_inference(ifos, dicc,targ_keys, priors, outdir, zenodo_file,folder):

    if dicc["inference"]["approximant"] == False:
        generated_waveform = waveform_generator(dicc,targ_keys)
        label = "DANSur"

    elif dicc["inference"]["approximant"] == True:
        generated_waveform = approximant_generator(dicc,targ_keys)
        label = "IMR"

    else:
        print('Approximant (in config.py) must either be True or False.')
        sys.exit(1)

    # ACTUALIZAR NUMPY dansur-22 1.0.0 requires numpy>=2.3.4, but you have numpy 1.26.4 which is incompatible.
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
    


    result=bilby.run_sampler(
                likelihood=likelihood,
                priors=priors,
                sampler="nessai",
                use_ratio=False,
                flow_class="gwflowproposal",
                npoints=dicc["inference"]["npoints"],
                resume=dicc["inference"]["resume"],
                outdir=outdir,
                label=label,
                npool=1,
                stopping=dicc["inference"]["stopping"]
            )

    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    plt.title("Reconstructed waveform", fontsize=14)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, 'reconstructed_waveform.png' ))
    return result



