import bilby
from bilby.gw.waveform_generator import LALCBCWaveformGenerator
def approximant_generator(dicc,targ_keys):
        waveform_arguments = dict(
            waveform_approximant=dicc["surrogate_model"], 
        )
        waveform_generator =bilby.gw.waveform_generator.WaveformGenerator(
            duration=dicc["duration"],
            sampling_frequency=dicc["sampling-frequency"],
            frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=waveform_arguments,
            start_time = dicc["start_time"] 
        )
        return waveform_generator

