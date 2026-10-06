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
sxs_simulation = "BBH:2157" 


# ============================================================
# 3. Inference configuration-
# ============================================================

resume = True          

# Supported models:            Link to respective papers and/or githubs):
#   non-precessing:
#     - DANSur                https://arxiv.org/abs/2412.06946 ; https://github.com/osvaldogramaxo/DANSur_22/
#     - IMRPhenomXHM
#     - NRHybSur3dq8          https://github.com/sxs-collaboration/gwsurrogate/blob/master/tutorial/website/NRHybSur2dq15.ipynb
#   precessing:
#     - NRSur7dq4             https://arxiv.org/abs/1905.09300

surrogate_model = 'dansur'

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

#Lista de eventos significativos:

# SXS:BBH:0110 ,  SXS:BBH:0237. non-precessing and asymetric simulations , https://arxiv.org/pdf/2012.11923 
# SXS:BBH:0926,  For this system, thengeneral relativistic phenomenon of spin-induced orbital precession42 is significant,  
#  This simulation was chosen since the majority
#of GW models obtain biased results, and disagree on the inferred binary parameters          https://eprints.soton.ac.uk/id/eprint/495815/1/2409.19404v1.pdf
#  SXS:BBH:1156 was chosen since it has largely asymmetric mass  https://eprints.soton.ac.uk/id/eprint/495815/1/2409.19404v1.pdf
#   SXS:BBH:1355 Igual masa, sin spin,  We consider a set of 3 public eccentric NR simulations
#SXS:BBH:1355, SXS:BBH:1359, and SXS:BBH:1363, which
#correspond to equal-mass, nonspinning BH binaries with initial eccentricities of 0.07, 0.13, and 0.25, respectively  https://arxiv.org/pdf/2412.12823
# These simulations have been assigned
#the identifiers SXS:BBH:1419 - SXS:BBH:1509,  https://arxiv.org/pdf/1812.07865
#


#We compare the recovered posteriors for 𝑚 1 , 𝑚 2 , 𝑞 , 𝜒 e f f m 1 ​ ,m 2 ​ ,q,χ eff ​ , the 90% credible intervals, and the log-evidence ln ⁡ 𝑍 lnZ for each model.”
#
