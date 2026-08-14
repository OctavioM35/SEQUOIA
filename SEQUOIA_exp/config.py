# config.py

results = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/gwtc_4_ruido_exp/pruebas_dansur"
data_folder = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/gwtc_4_ruido_exp/gwtx_4_ruido"

#Infer a real event
run_particular_event = True             #True for running one event, False for running all events on data_folder directory
particular_event = 'GW231001_140220'
resume = False                          #True for resuming one event if stopped previously, False for running the event from the begining
approximant = True                      #True: IMR estimates waveform , False: DANSur estimates the waveform
skip = False
npoints = 1000
stopping = 1


#Generation and inference of a synthetic singal:
output_synthetic = '/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/sintetico'
generate_synthetic_signal = False
injection_parameters = {
                "chirp_mass": 30,
                "mass_ratio": 0.5,
                "chi_1": 0.0,
                "chi_2":0.0,
                "tilt_1": 0.0,
                "tilt_2": 0.0,
                "phi_12": 0.0,
                "phi_jl": 0.0,
                "luminosity_distance": 4000,
                "theta_jn":0,
                "psi": 0,
                "phase": 0.0,
                "geocent_time": 0,
                "ra":0,
                "dec": 0,
            }


automatic = True                        #True for using all Zenodo priors, duration and sampling frecuencies. False for customs
duration = 8
sampling_frequency = 1024               #Hz

dicc = {                                # Modify for custom configuration only if automatic == False
    # Data configuration
    "duration": 9,
    "sampling_frequency": 1024,

    # Frequency configuration
    "frequencies": {
        "h1": {
            "minimum": 20,
            "maximum": 480
        },
        "l1": {
            "minimum": 20,
            "maximum": 489
        },
        "v1": {
            "minimum": 20,
            "maximum": 480
        },
        
    },

    # Waveform configuration
    "waveform": {
        "reference-frequency": 0,
        "minimum-frequency": None,
        "pn-spin-order": 0,
        "pn-phase-order": 0,
        "pn-tidal-order": 0,
        "pn-amplitude-order": 0,
        "mode-array": None,
        "catch-waveform-errors": None,
        "f_ref":0,
        "f_low":0
    },

    # Event configuration
    "gps_time": 0,

    # Marginalization configuration
    "marginalization": {
        "distance": False,
        "phase": True,
        "time": True,
        "jitter-time": False,
        "calibration": None,
    },

    # Inference configuration
    "inference": {
        "resume": resume,
        "npoints": npoints,
        "stopping": stopping,
        "approximant": approximant,
        "skip": skip,
        "sampling-seed":0,
    },

    # General configuration
    "run_particular_event": run_particular_event,
    "automatic": automatic,
    "generate_synthetic_signal": generate_synthetic_signal,
    "injection_parameters": injection_parameters
    
}


# MIRAR EN DETALLE ESTOS EVENTOS (SON LOS RUIDOSOS)
# ['GW230628_231200', 'GW230819_171910', 'GW231206_233901', 'GW230825_041334', 'GW240109_050431', 'GW231029_111508', 'GW240107_013215', 'GW230814_230901', 'GW230820_212515', 'GW230601_224134', 'GW230803_033412', 'GW231001_140220', 'GW231004_232346', 'GW231113_122623', 'GW230914_111401', 'GW230814_061920']


