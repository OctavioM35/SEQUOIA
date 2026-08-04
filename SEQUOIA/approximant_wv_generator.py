import bilby
from bilby.gw.waveform_generator import LALCBCWaveformGenerator
def approximant_generator(dicc,targ_keys):

    waveform_arguments = dict(
        waveform_approximant="IMRPhenomXO4a",
        minimum_frequency=int(dicc["minimum-frequency"][" waveform"]),
        pn_spin_order=dicc["pn-spin-order"],
        pn_phase_order=dicc["pn-phase-order"],
        pn_tidal_order=dicc["pn-tidal-order"],
        pn_amplitude_order=dicc["pn-amplitude-order"],
        mode_array=dicc["mode-array"],
        catch_waveform_errors=dicc["catch-waveform-errors"],
        reference_frequency=dicc["reference-frequency"],
    )

    waveform_generator = LALCBCWaveformGenerator(
        duration=dicc["duration"],
        sampling_frequency=dicc["sampling-frequency"],
        frequency_domain_source_model=(bilby.gw.source.lal_binary_black_hole),
        parameter_conversion=(bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters),
        waveform_arguments=waveform_arguments,
    )
    waveform_generator.source_parameter_keys = targ_keys

    return waveform_generator

