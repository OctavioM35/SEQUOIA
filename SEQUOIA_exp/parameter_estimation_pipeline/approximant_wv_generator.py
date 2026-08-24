import bilby
from bilby.gw.waveform_generator import LALCBCWaveformGenerator
def approximant_generator(dicc,targ_keys):
        waveform_arguments = dict(
            waveform_approximant="IMRPhenomXO4a", #Versión con precesión: IMRPhenomXO4a  #Versión sin precesión IMRPhenomXHM
            # minimum_frequency=int(dicc["waveform"]["minimum-frequency"]),
            # pn_spin_order=dicc["waveform"]["pn-spin-order"],
            # pn_phase_order=dicc["waveform"]["pn-phase-order"],
            # pn_tidal_order=dicc["waveform"]["pn-tidal-order"],
            # pn_amplitude_order=dicc["waveform"]["pn-amplitude-order"],
            # mode_array=dicc["waveform"]["mode-array"],
            # catch_waveform_errors=dicc["waveform"]["catch-waveform-errors"],
            # reference_frequency=dicc["waveform"]["reference-frequency"],
        )
        waveform_generator =bilby.gw.waveform_generator.WaveformGenerator(
            duration=dicc["duration"],
            sampling_frequency=dicc["sampling-frequency"],
            frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=waveform_arguments,
            start_time = dicc["start_time"] 
        )
        # waveform_generator.source_parameter_keys = set(targ_keys)

        return waveform_generator

