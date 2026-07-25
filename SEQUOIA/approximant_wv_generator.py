import bilby

def approximant_generator(duration,sampling_frequency):    
    waveform_arguments = dict(
        waveform_approximant="IMRPhenomPv2",  # aproximant nativo, es el modelo que generará el waveform
        minimum_frequency=20.0  #EL likelihood no utilizará parte de la señal con w<minimum_frequency
    )

    waveform_generator = bilby.gw.waveform_generator.LALCBCWaveformGenerator(  
        duration=duration,                                                      
        sampling_frequency=sampling_frequency,                                 
        frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,    
        # time_domain_source_model=bilby.gw.source.lal_binary_black_hole,
        parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters, 
        waveform_arguments=waveform_arguments                                  

    )
    return waveform_generator