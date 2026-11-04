import ast
import h5py
import numpy as np
from config import (
    # surrogate_model,
    # run_particular_event,
    resume,
    npoints,
    automatic,
    stopping,
    skip,
)



# ---------------------------------------------------------
# 5. Load LVK configuration
# ---------------------------------------------------------
def print_dict(d, indent=0):
    """Print a nested dictionary as an indented list of "- key: value" lines."""

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
    Read a single-element HDF5 dataset and return a Python value.

    Bytes are decoded to strings, strings are evaluated with
    ast.literal_eval when possible (otherwise returned as text), and numpy
    scalars are converted to native Python types.
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
    """
    Create the default configuration dictionary of the inference.

    The dictionary groups the settings used along the pipeline: data
    (duration and sampling frequency), minimum and maximum frequencies per
    detector, waveform settings (reference, low, start and final
    frequencies, PN orders, mode array), event times (GPS, start,
    geocentric and trigger time, Tukey roll-off, post-trigger), the
    marginalizations, the sampler options (resume, npoints, stopping, skip,
    seed) and general options (automatic mode, precession, surrogate model
    and label). Frequencies set to 0 mean "not specified". The values are
    later overwritten with those read from the LVK Zenodo file.

    INPUTS:
    -------

    -None. The sampler options (resume, npoints, stopping, skip) and the
    automatic option are taken from config.py.

    OUTPUTS:
    --------

    -dicc: New dictionary with the default configuration.
    """

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
        "precession_model": True,
        "surrogate_model": None,
        "label": 'to be determined',

    }

    return dicc


def rewrite(s):
    """
    Convert a string representing a dictionary into a Python dictionary.

    The function removes the braces, splits the string into "key: value"
    pairs and evaluates each value with ast.literal_eval. If a value cannot
    be evaluated, it is kept as a string. It is used to read the
    per-detector frequencies stored in the Zenodo configuration, for
    example "{'H1': 20, 'L1': 20, 'V1': 20}".

    INPUTS:
    -------

    -s: String with the dictionary stored in the Zenodo file.

    OUTPUTS:
    --------

    -Dictionary with the keys (as strings, quotes included) and the
    evaluated values.

    -The original string s if a pair does not contain a ":" separator.
    """

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
    """
    Load the configuration used by the LVK for an event from its Zenodo
    file.

    The function starts from the default configuration and overwrites it
    with the values stored in the first analysis of the file. From the
    "config_file" group it reads the duration, the minimum, maximum and
    reference frequencies (per detector when given as a dictionary), the
    PN orders, the mode array, the marginalizations, the roll-off,
    the post-trigger time, the spline calibration nodes, the sampling seed
    and the geocentric, trigger, GPS and start times. From the "meta_data"
    group it reads the trigger, GPS and geocentric times, the start time
    and the f_ref, f_low, f_start and f_final frequencies, which take
    precedence over the previous ones. The final configuration is printed.

    INPUTS:
    -------

    -zenodo_file: Path to the LVK Zenodo .h5 file of the event.

    OUTPUTS:
    --------

    -dicc: Configuration dictionary with the values of the LVK analysis.
    Parameters not present in the file keep their default value.
    """

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
