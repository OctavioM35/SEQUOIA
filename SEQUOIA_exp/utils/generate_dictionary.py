import ast
import h5py

from config import (
    run_particular_event,
    resume,
    approximant,
    npoints,
    automatic,
    generate_synthetic_signal,
    injection_parameters,
    stopping,
    skip,
)

def get_default_configuration():

    dicc = {

        # =====================================================
        # Data configuration
        # =====================================================

        "duration": 8,
        "sampling-frequency": 1024,

        # =====================================================
        # Frequency configuration
        # =====================================================

        "frequencies": {

            "h1": {
                "minimum": 0,
                "maximum": 0,
            },

            "l1": {
                "minimum": 0,
                "maximum": 0,
            },

            "v1": {
                "minimum": 0,
                "maximum": 0,
            },
        },

        # =====================================================
        # Waveform configuration
        # =====================================================

        "waveform": {
            "reference-frequency": 0,
            "minimum-frequency": 0,
            "pn-spin-order": 0,
            "pn-phase-order": 0,
            "pn-tidal-order": 0,
            "pn-amplitude-order": 0,
            "mode-array": None,
            "catch-waveform-errors": None,
            "f_ref":0,
            "f_low":0,
            "f_start":0,
            "f_final": 0,
        },

        # =====================================================
        # Event configuration
        # =====================================================

        "gps_time": 0,

        # =====================================================
        # Marginalization configuration
        # =====================================================

        "marginalization": {

            "distance": None,
            "phase": None,
            "time": None,
            "jitter": None,
            "calibration": None,
        },

        # =====================================================
        # Inference configuration
        # =====================================================

        "inference": {

            "resume": resume,
            "npoints": npoints,
            "stopping": stopping,
            "approximant": approximant,
            "skip": skip,
        },

        # =====================================================
        # General configuration
        # =====================================================

        "automatic": automatic,

        "generate_synthetic_signal": generate_synthetic_signal,

        # =====================================================
        # Injection configuration
        # =====================================================

        "injection_parameters": (
            injection_parameters
            if generate_synthetic_signal
            else None
        ),
    }

    return dicc


def rewrite(s):
    d={}
    config = s.strip('{}  ,').split(',')
    for i in config:
            param=i.split(':')
            try:
                d[param[0]] = ast.literal_eval(param[1])
            except (ValueError, SyntaxError):
                d[param[0]] = param[1]
            except IndexError:
                return s    
    return d

def load_zenodo_configuration(zenodo_file):
    dicc = get_default_configuration()

    with h5py.File(zenodo_file, "r") as f:
            key0 = list(f.keys())[0]
            print('=' *50)
            print('Loading Zenodo configuration...')
            print(' ')
            config = f[key0]["config_file"]["config"]
            envelope = f[key0]["calibration_envelope"]
            meta=f[key0]["meta_data"]["meta_data"]
            # for data in (config,meta):
            for parameter, ds in config.items():

                    if parameter in dicc:
                        s = ds[0].decode("utf-8").strip()
                        dicc[parameter] = ast.literal_eval(s)


                    if "marginalization" in parameter or "jitter" in parameter:
                        s = ds[0].decode("utf-8").strip()
                        try:
                            param=parameter.split("-")
                            dicc["marginalization"][param[0]] = ast.literal_eval(s)
                        except ValueError:
                            pass


                    if "frequency" in parameter:
                        s = ds[0].decode("utf-8").strip()
                        param = parameter.split("-")

                        if param[0] =='reference':
                            dicc["waveform"][parameter] = ast.literal_eval(s)

                        elif param[0] == 'sampling':
                            pass

                        elif s.startswith('{') and s.endswith('}'):
                            frequency_dictionary = rewrite(s)

                            for detector_frequency in  frequency_dictionary:

                                if param[0] =='maximum':
                                    if detector_frequency =='H1':
                                        dicc["frequencies"]["h1"]["maximum"] = frequency_dictionary[detector_frequency]
                                    if detector_frequency =='L1':
                                        dicc["frequencies"]["l1"]["maximum"] = frequency_dictionary[detector_frequency]
                                    if detector_frequency =='V1':
                                        dicc["frequencies"]["v1"]["maximum"] = frequency_dictionary[detector_frequency]

                                          
                                if param[0] =='minimum':
                                    if detector_frequency =='H1':
                                        dicc["frequencies"]["h1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if detector_frequency =='L1':
                                        dicc["frequencies"]["l1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if detector_frequency =='V1':
                                        dicc["frequencies"]["v1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if detector_frequency =='waveform':
                                        dicc["waveform"]["minimum-frequency"] = frequency_dictionary[detector_frequency]

                    if 'pn' in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc['waveform'][parameter] = ast.literal_eval(s)

                    if parameter =='mode-array' or parameter =='catch-waveform-errors':
                        s = ds[0].decode("utf-8").strip()
                        dicc['waveform'][parameter]= ast.literal_eval(s)

                    if parameter =='trigger-time' or parameter =='t_gps' or parameter =='geocent_time':
                        s = ds[0].decode("utf-8").strip()
                        dicc["gps_time"] = ast.literal_eval(s)     
                    if parameter =='f_ref' or parameter =='f_low' or parameter =='f_start' or parameter =='f_final': 
                        s = ds[0].decode("utf-8").strip()                                            
                        dicc["waveform"][parameter] = ast.literal_eval(s)     

            for parameter in meta.items():
                    if parameter =='trigger-time' or parameter =='t_gps' or parameter =='geocent_time':
                        dicc["gps_time"] = ast.literal_eval(s)     
                    if parameter =='f_ref' or parameter =='f_low' or parameter =='f_start' or parameter =='f_final': 
                        dicc["waveform"][parameter] = ast.literal_eval(s)     
            print('Loading completed.')

            return dicc
