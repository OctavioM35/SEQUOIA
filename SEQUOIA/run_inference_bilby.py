import bilby
from bilby.gw.likelihood import GravitationalWaveTransient
from waveform_generator import  waveform_generator
from matplotlib import pyplot as plt
import os
from gwpy.timeseries import TimeSeries


def run_inference(ifos, duration, sampling_frequency,targ_keys, priors, npoints, outdir,resume):
    label = "NNSur"
    generated_waveform = waveform_generator(duration, sampling_frequency,targ_keys)

    likelihood = GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=True,
        time_marginalization=True,
        jitter_time=False
    )


    result=bilby.run_sampler(
            likelihood=likelihood,
            priors=priors,
            sampler="nessai",
            use_ratio=False,
            flow_class="gwflowproposal",
            npoints=npoints,
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
    plt.savefig(os.path.join(outdir, 'reconstructed_waveform.png'))
    return result