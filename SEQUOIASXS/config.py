# config.py

# ============================================================
# 1. Paths
# ============================================================

results = (
    "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/simulations/NRH"
)

# ============================================================
# 2. Event selection
# ============================================================

#True -> runs only particular_event
#False -> runs all events on data folder
run_particular_event = False     

# Event to analyse when run_particular_event = True
sxs_simulation = "BBH:0926" 


# ============================================================
# 3. Inference configuration-
# ============================================================

resume = False          

# Supported models:            Link to respective papers and/or githubs):
#   non-precessing:
#     - DANSur                https://arxiv.org/abs/2412.06946 ; https://github.com/osvaldogramaxo/DANSur_22/
#     - IMRPhenomXHM
#     - NRHybSur3dq8          https://github.com/sxs-collaboration/gwsurrogate/blob/master/tutorial/website/NRHybSur2dq15.ipynb
#   precessing:
#     - NRSur7dq4             https://arxiv.org/abs/1905.09300

surrogate_model = 'nrh'

#False -> reruns already finished events
#True -> skips processed events
skip = True

npoints = 1000

stopping = 0.1

sampling_seed = 0

# ============================================================
# 6. dicc initialization
# ============================================================

dicc = {

    # --------------------------------------------------------
    # Duration and sampling frequency configuration
    # --------------------------------------------------------

    "duration": 3,

    "sampling-frequency": 4096,

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
    "surrogate_model": surrogate_model,
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
    },
        # =====================================================
        # Inference configuration
        # =====================================================

        "inference": {

            "resume": resume,
            "npoints": npoints,
            "stopping": stopping,
            "skip": skip,
            "sampling-seed": sampling_seed,

        },

        "signal": {
        "timeseries_H1": 0,
        "h_plus_uniform": 0,
        "h_cross_uniform": 0,
        "t_uniform": 0,
         },
        "signal_sxs": None,
        "generated_waveform": None,
        "inyection" : None
    }


