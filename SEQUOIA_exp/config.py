# config.py

# ============================================================
# Paths
# ============================================================

results = (
    "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/resultados_paper"
)

data_folder = (
    "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/gwtc_9s"
)

# ============================================================
# Event selection
# ============================================================
#True -> runs only particular_event
#False -> runs all events on data folder
run_particular_event = False

# Event to analyse when run_particular_event = True
particular_event = "GW191204_110529"


# ============================================================
# Inference configuration
# ============================================================

resume = True

# False -> use DANSur waveform
# True  -> use IMR waveform
approximant = False

skip = True

npoints = 1000

stopping = 0.1

sampling_seed = 0

# ============================================================
# Analysis mode
# ============================================================

# True  -> use Zenodo configuration
# False -> use custom configuration below
automatic = True


# ============================================================
# Synthetic signal configuration
# ============================================================

generate_synthetic_signal = False

# Used only when generate_synthetic_signal = True
injection_parameters = {

    "chirp_mass": 30,
    "mass_ratio": 0.5,

    "chi_1": 0.0,
    "chi_2": 0.0,

    "tilt_1": 0.0,
    "tilt_2": 0.0,

    "phi_12": 0.0,
    "phi_jl": 0.0,

    "luminosity_distance": 4000,

    "theta_jn": 0.0,

    "psi": 0.0,
    "phase": 0.0,

    "geocent_time": 0,

    "ra": 0.0,
    "dec": 0.0,
}


# ============================================================
# Custom configuration
#
# Used only when automatic = False
# ============================================================

custom_config = {

    # --------------------------------------------------------
    # Data
    # --------------------------------------------------------

    "duration": 9,

    "sampling_frequency": 1024,


    # --------------------------------------------------------
    # Frequencies
    # --------------------------------------------------------

    "frequencies": {

        "h1": {
            "minimum": 20,
            "maximum": 480,
        },

        "l1": {
            "minimum": 20,
            "maximum": 489,
        },

        "v1": {
            "minimum": 20,
            "maximum": 480,
        },
    },


    # --------------------------------------------------------
    # Waveform
    # --------------------------------------------------------

    "waveform": {

        "reference-frequency": 0,

        "minimum-frequency": None,

        "pn-spin-order": 0,
        "pn-phase-order": 0,
        "pn-tidal-order": 0,
        "pn-amplitude-order": 0,

        "mode-array": None,

        "catch-waveform-errors": None,
    },


    # --------------------------------------------------------
    # Event
    # --------------------------------------------------------

    "gps_time": 0,


    # --------------------------------------------------------
    # Marginalization
    # --------------------------------------------------------

    "marginalization": {

        "distance": False,
        "phase": True,
        "time": True,
        "jitter-time": False,
        "calibration": None,
    },
}

# MIRAR EN DETALLE ESTOS EVENTOS (SON LOS RUIDOSOS)
# ['GW230628_231200', 'GW230819_171910', 'GW231206_233901', 'GW230825_041334', 'GW240109_050431', 'GW231029_111508', 'GW240107_013215', 'GW230814_230901', 'GW230820_212515', 'GW230601_224134', 'GW230803_033412', 'GW231001_140220', 'GW231004_232346', 'GW231113_122623', 'GW230914_111401', 'GW230814_061920']


