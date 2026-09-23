# config.py

# ============================================================
# 1. Paths
# ============================================================

results = (
    "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master"
)

data_folder = (
    "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/gwtc_selected"
)

# ============================================================
# 2. Event selection
# ============================================================

#True -> runs only particular_event
#False -> runs all events on data folder
run_particular_event = True     

# Event to analyse when run_particular_event = True
particular_event = "GW200216_220804" 
# ============================================================
# 3. Inference configuration
# ============================================================

resume =True

# Supported models:            Link to respective papers and/or githubs:
#   non-precessing:
#     - DANSur                https://arxiv.org/abs/2412.06946 ; https://github.com/osvaldogramaxo/DANSur_22/
#     - IMRPhenomXHM          https://arxiv.org/abs/2001.10914
#     - NRHybSur3dq8          https://github.com/sxs-collaboration/gwsurrogate/blob/master/tutorial/website/NRHybSur2dq15.ipynb
#   precessing:
#     - IMRPhenomXO4a         https://arxiv.org/abs/2312.10025
#     - NRSur7dq4             https://arxiv.org/abs/1905.09300

surrogate_model = 'dansur'

#False -> reruns already finished events
#True -> skips processed events
skip = False

npoints =500

stopping = 1

sampling_seed = 0

# ============================================================
# 4. Analysis mode
# ============================================================

# True  -> use Zenodo configuration  (duration, frequencies, marginalization, priors, times)
# False -> use custom configuration (customize in section 6)
automatic = True


# ============================================================
# 5. Synthetic signal configuration

# It should be noted that automatic == True effects will 
# still apply to synthetic signals.
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
# 6. Custom configuration
#
# Used only when automatic = False
# ============================================================

custom_config = {

    # --------------------------------------------------------
    # Duration and sampling frequency configuration
    # --------------------------------------------------------

    "duration": 9,

    "sampling_frequency": 1024,


    # --------------------------------------------------------
    # Interferometer maximum and minimum frequencies configuration
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
    # Waveform parameter configuration
    # --------------------------------------------------------

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


    # --------------------------------------------------------
    # Event times configuration
    # --------------------------------------------------------

    "gps_time": 0,
    "start_time": 0,
    "geocent_time": 0,
    "trigger_time":  0,
    "roll_off": 0.2,

    # --------------------------------------------------------
    # Marginalization configuration
    # --------------------------------------------------------

    "marginalization": {

        "distance": False,
        "phase": True,
        "time": True,
        "jitter-time": False,
        "calibration": None,
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
        # =====================================================
        # Injection configuration
        # =====================================================

        "injection_parameters": (
            injection_parameters
            if generate_synthetic_signal
            else None
        ),
    }
