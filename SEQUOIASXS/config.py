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
run_particular_event = True     

# Event to analyse when run_particular_event = True
sxs_simulation = "BBH:0198" 

# ============================================================
# 3. Inference configuration-
# ============================================================

resume = False           

# Supported models:            Link to respective papers and/or githubs):
#     - DANSur                https://arxiv.org/abs/2412.06946 ; https://github.com/osvaldogramaxo/DANSur_22/
#     - NRHybSur3dq8          https://github.com/sxs-collaboration/gwsurrogate/blob/master/tutorial/website/NRHybSur2dq15.ipynb
#     - NRSur7dq4             https://arxiv.org/abs/1905.09300

surrogate_model = 'dansur'

#False -> reruns already finished events
#True -> skips processed events
skip = False

npoints = 500

stopping = 4

sampling_seed = 0

# ============================================================
# 4. dicc initialization
# ============================================================
dicc = {
    "geocent_time": 'to be determined',

    # --------------------------------------------------------
    # Duration and sampling frequency configuration
    # --------------------------------------------------------

    "duration": 'to be determined',

    "sampling-frequency": 4096,

    # --------------------------------------------------------
    # Waveform parameter configuration
    # --------------------------------------------------------

    "waveform": {

        "reference-frequency": 'to be determined',
        "minimum-frequency": 'to be determined',

    },


    # --------------------------------------------------------
    # Event times configuration
    # --------------------------------------------------------
    "surrogate_model": surrogate_model,
    "geocent_time": 'to be determined',
    "trigger_time":  'to be determined',
    "roll_off": 0,

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
        "timeseries_H1": 'to be determined',
        "h_plus_uniform": 'to be determined',
        "h_cross_uniform": 'to be determined',
        "t_uniform": 'to be determined',
         },
        "signal_sxs": None,
        "generated_waveform": None,
        "inyection" : None,
        "phase": None
    }

#Lista de eventos significativos:

# SXS:BBH:0110 ,  SXS:BBH:0237. non-precessing and asymetric simulations , https://arxiv.org/pdf/2012.11923 
# SXS:BBH:0926,  For this system, thengeneral relativistic phenomenon of spin-induced orbital precession42 is significant,  
#  This simulation was chosen since the majority
#of GW models obtain biased results, and disagree on the inferred binary parameters          https://eprints.soton.ac.uk/id/eprint/495815/1/2409.19404v1.pdf
#  SXS:BBH:1156 was chosen since it has largely asymmetric mass  https://eprints.soton.ac.uk/id/eprint/495815/1/2409.19404v1.pdf
#   SXS:BBH:1355 Igual masa, sin spin,  We consider a set of 3 public eccentric NR simulations
#, SXS:BBH:1359, and SXS:BBH:1363, which
#correspond to equal-mass, nonspinning BH binaries with initial eccentricities of 0.07, 0.13, and 0.25, respectively  https://arxiv.org/pdf/2412.12823
# These simulations have been assigned
#the identifiers SXS:BBH:1419 - SXS:BBH:1509,  https://arxiv.org/pdf/1812.07865
#

