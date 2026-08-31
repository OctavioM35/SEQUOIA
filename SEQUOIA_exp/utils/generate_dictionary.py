import ast
import h5py
import numpy as np
from config import (
    surrogate_model,
    run_particular_event,
    resume,
    npoints,
    automatic,
    generate_synthetic_signal,
    injection_parameters,
    stopping,
    skip,
)



# ---------------------------------------------------------
# 5. Load LVK configuration
# ---------------------------------------------------------
def print_dict(d, indent=0):
    for key, value in d.items():

        spaces = " " * indent

        if isinstance(value, dict):
            print(f"{spaces}- {key}:")
            print_dict(value, indent + 4)

        elif isinstance(value, list):
            print(f"{spaces}- {key}:")
            for item in value:
                if isinstance(item, dict):
                    print_dict(item, indent + 4)
                else:
                    print(f"{spaces}    - {item}")

        else:
            print(f"{spaces}- {key}: {value}")

def read_hdf5_value(value):
    """
    Lee un dataset HDF5 de un solo elemento y devuelve un valor Python.
    Funciona para bytes, strings, ints, floats y tipos numpy.
    """

    raw = value[0]

    # bytes / numpy.bytes_
    if isinstance(raw, (bytes, np.bytes_)):
        raw = raw.decode("utf-8").strip()

    # String
    if isinstance(raw, str):
        try:
            return ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            return raw

    # Tipos numpy -> tipos Python
    if isinstance(raw, np.generic):
        return raw.item()

    # int, float, bool, etc.
    return raw


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
        "start_time": 0,
        "geocent_time": 0,
        "trigger_time":  0,
        "roll_off": 0.2,
        "post_trigger": 0,
        "time_reference": 'gecentric',
        "reference_frame": 'H1L1',
        "spline-calibration-nodes":20,
        # =====================================================
        # Marginalization configuration
        # =====================================================

        "marginalization": {

            "distance": None,
            "phase": None,
            "time": None,
            "jitter": None,
            "calibration": True,
        },

        # =====================================================
        # Inference configuration
        # =====================================================

        "inference": {

            "resume": resume,
            "npoints": npoints,
            "stopping": stopping,
            "skip": skip,
            "sampling-seed": None,

        },

        # =====================================================
        # General configuration
        # =====================================================

        "automatic": automatic,
        "generate_synthetic_signal": generate_synthetic_signal,
        "precession_model": True,
        "surrogate_model": surrogate_model,
        "label": 0,
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

                    if parameter in dicc and parameter != 'label':
                        s = ds[0].decode("utf-8").strip()
                        dicc[parameter] = ast.literal_eval(s)


                    if "marginalization" in parameter or "jitter" in parameter:
                        s = ds[0].decode("utf-8").strip()
                        try:
                            param=parameter.split("-")
                            dicc["marginalization"][param[0]] = ast.literal_eval(s)
                        except ValueError:
                            pass
                        except SyntaxError:
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

                                if param[0] == 'maximum':
                                    if 'H1' in detector_frequency:
                                        dicc["frequencies"]["h1"]["maximum"] = frequency_dictionary[detector_frequency]
                                    if 'L1' in detector_frequency:
                                        dicc["frequencies"]["l1"]["maximum"] = frequency_dictionary[detector_frequency]
                                    if 'V1' in detector_frequency:
                                        dicc["frequencies"]["v1"]["maximum"] = frequency_dictionary[detector_frequency]


                                if param[0] == 'minimum':
                                    if 'H1' in detector_frequency:
                                        dicc["frequencies"]["h1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if 'L1' in detector_frequency:
                                        dicc["frequencies"]["l1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if 'V1' in detector_frequency:
                                        dicc["frequencies"]["v1"]["minimum"] = frequency_dictionary[detector_frequency]
                                    if 'waveform' in detector_frequency:
                                        dicc["waveform"]["minimum-frequency"] = frequency_dictionary[detector_frequency]


                    if 'pn' in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc['waveform'][parameter] = ast.literal_eval(s)

                    if 'geocent-time' in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc["geocent_time"] = ast.literal_eval(s)

                    if parameter =='mode-array' or parameter =='catch-waveform-errors':
                        s = ds[0].decode("utf-8").strip()
                        dicc['waveform'][parameter]= ast.literal_eval(s)

                    if  parameter =='t_gps':
                        s = ds[0].decode("utf-8").strip()
                        dicc["gps_time"] = ast.literal_eval(s)     

                    if parameter =='f_ref' or parameter =='f_low' or parameter =='f_start' or parameter =='f_final': 
                        s = ds[0].decode("utf-8").strip()                             
                        dicc["waveform"][parameter] = ast.literal_eval(s)     

                    if parameter =='sampling-seed':
                        s = ds[0].decode("utf-8").strip()                                            
                        dicc["inference"][parameter] = ast.literal_eval(s)     

                    if 'start_time' == parameter:
                        s = ds[0].decode("utf-8").strip()                                            
                        dicc["start_time"] = ast.literal_eval(s)     

                    if parameter == "trigger-time":
                        s = ds[0].decode("utf-8").strip()
                        dicc["trigger_time"] =  ast.literal_eval(s)

                    if parameter == "geocent_time":
                        s = ds[0].decode("utf-8").strip()
                        dicc["geocent_time"] =  ast.literal_eval(s)

                    if "roll-off" in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc["roll_off"] =  ast.literal_eval(s)

                    if "post-trigger" in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc["post_trigger"] =  ast.literal_eval(s)
                    if 'spline-calibration-nodes' in parameter:
                        s = ds[0].decode("utf-8").strip()
                        dicc["spline-calibration-nodes"] =  ast.literal_eval(s)

                    # if "reference-frame" in parameter:
                    #     s = ds[0].decode("utf-8").strip()
                    #     dicc["reference_frame"] =  ast.literal_eval(s)

                    # if "time-reference" in parameter:
                    #   s = ds[0].decode("utf-8").strip()
                    #   dicc["time_reference"] =  ast.literal_eval(s)

            for parameter, value in meta.items():
                if parameter == "trigger-time":
                    dicc["trigger_time"] = read_hdf5_value(value)

                if parameter == "t_gps":
                    dicc["t_gps"] = read_hdf5_value(value)

                if parameter == "geocent_time":
                    dicc["geocent_time"] = read_hdf5_value(value)

                elif parameter in ["f_ref", "f_low", "f_start", "f_final"]:
                    dicc["waveform"][parameter] = read_hdf5_value(value)

    
                elif parameter == "start_time":

                    dicc["start_time"] = read_hdf5_value(value)


            print("Loading completed. Parameters from LVK loaded:\n")
            print_dict(dicc)
            return dicc
