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

sxs_simulation = [ 'SXS:BBH:0032', 'SXS:BBH:0045', 'SXS:BBH:1412', 'SXS:BBH:0007']
# BBH:0198 , BBH:0237

# ============================================================
# 3. Inference configuration-
# ============================================================

resume = True           

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
# Nueva lisrta de eventos significativos: 
#  SXS:BBH:0237, SXS:BBH:0198 ,  SXS:BBH:1942, SXS:BBH:1485, SXS:BBH:3630
# Total simulaciones seleccionadas: 212
#                  q_norm        chi1_z        chi2_z  reference_eccentricity
# SXS ID                                                                     
# SXS:BBH:0007 0.6667  0.0000  0.0000 0.000434  non_spinning
# SXS:BBH:0031 0.3332  0.4995  0.0000 0.000050 spinning_asym
# SXS:BBH:0036 0.3332 -0.4995  0.0000 0.000533 spinning_asym
# SXS:BBH:0040 0.3332 -0.4995  0.0000 0.000911 spinning_asym
# SXS:BBH:0045 0.3334  0.4995 -0.4994 0.000662 spinning_asym
# SXS:BBH:0174 0.3332  0.4995  0.0000 0.000295 spinning_asym
# SXS:BBH:0191 0.3988 -0.0000  0.0001 0.000758  non_spinning
# SXS:BBH:0194 0.6590  0.0000 -0.0001 0.000774  non_spinning
# SXS:BBH:0198 0.8316 -0.0001  0.0001 0.000202  non_spinning
# SXS:BBH:0200 0.3056 -0.0001 -0.0000 0.000394  non_spinning
# SXS:BBH:0201 0.4317  0.0001 -0.0000 0.000194  non_spinning
# SXS:BBH:0236 0.5000 -0.6000 -0.0000 0.000095 spinning_asym
# SXS:BBH:0244 0.5000 -0.0001 -0.5999 0.000034 spinning_asym
# SXS:BBH:0247 0.5000 -0.0000  0.6000 0.000132 spinning_asym
# SXS:BBH:0254 0.5000  0.6000 -0.5999 0.000052 spinning_asym
# SXS:BBH:0255 0.5000  0.6000 -0.0000 0.000038 spinning_asym
# SXS:BBH:0259 0.4000  0.0000  0.0000 0.000009  non_spinning
# SXS:BBH:0262 0.3333 -0.6000 -0.0000 0.000203 spinning_asym
# SXS:BBH:0263 0.3334 -0.6000  0.5998 0.000205 spinning_asym
# SXS:BBH:0266 0.3333 -0.6000  0.3999 0.000170 spinning_asym
# SXS:BBH:0269 0.3333 -0.4000  0.6000 0.000116 spinning_asym
# SXS:BBH:0275 0.3333 -0.0000 -0.5998 0.000106 spinning_asym
# SXS:BBH:0278 0.3334  0.0001  0.5999 0.000198 spinning_asym
# SXS:BBH:0284 0.3333  0.4000 -0.5998 0.000115 spinning_asym
# SXS:BBH:0287 0.3333  0.6000 -0.5998 0.000072 spinning_asym
# SXS:BBH:0288 0.3333  0.5999 -0.4000 0.000210 spinning_asym
# SXS:BBH:0289 0.3333  0.5999 -0.0000 0.000230 spinning_asym
# SXS:BBH:0294 0.2857  0.0000  0.0000 0.000042  non_spinning
# SXS:BBH:0310 0.8189  0.0001  0.0001 0.000558  non_spinning
# SXS:BBH:0317 0.3006  0.5226 -0.4483 0.000640 spinning_asym
# SXS:BBH:0331 0.5001 -0.0001 -0.7996 0.000154 spinning_asym
# SXS:BBH:0332 0.5001  0.0002  0.7995 0.000393 spinning_asym
# SXS:BBH:0354 0.5001 -0.2002  0.4000 0.000113 spinning_asym
# SXS:BBH:0355 0.5715 -0.6003  0.7997 0.000663 spinning_asym
# SXS:BBH:0371 0.5714 -0.0001 -0.4000 0.000837 spinning_asym
# SXS:BBH:0375 0.5001 -0.0001 -0.3999 0.000644 spinning_asym
# SXS:BBH:0382 0.5715 -0.0002  0.7999 0.000327 spinning_asym
# SXS:BBH:0388 0.5713  0.7996  0.4001 0.000517 spinning_asym
# SXS:BBH:0407 0.5001 -0.0001  0.4000 0.000724 spinning_asym
# SXS:BBH:0410 0.5001  0.6001 -0.0001 0.000036 spinning_asym
# SXS:BBH:0414 0.5715 -0.4002 -0.7997 0.000242 spinning_asym
# SXS:BBH:0423 0.5715  0.0001 -0.7996 0.000267 spinning_asym
# SXS:BBH:0435 0.5715  0.4001  0.7999 0.000106 spinning_asym
# SXS:BBH:0448 0.5000  0.4001 -0.4001 0.000124 spinning_asym
# SXS:BBH:0451 0.5715  0.0002  0.3998 0.000087 spinning_asym
# SXS:BBH:0461 0.5000 -0.6005  0.0001 0.000180 spinning_asym
# SXS:BBH:0473 0.5715  0.2001 -0.7996 0.000212 spinning_asym
# SXS:BBH:0488 0.5715  0.6001 -0.7996 0.000464 spinning_asym
# SXS:BBH:0513 0.5000  0.6000 -0.4002 0.000457 spinning_asym
# SXS:BBH:0550 0.5715 -0.2003  0.7999 0.000127 spinning_asym
# SXS:BBH:0552 0.5714  0.7999 -0.4000 0.000079 spinning_asym
# SXS:BBH:0554 0.5000  0.2001 -0.4001 0.000101 spinning_asym
# SXS:BBH:0566 0.5715  0.1999  0.7999 0.000398 spinning_asym
# SXS:BBH:0593 0.6667  0.0000  0.0001 0.000058  non_spinning
# SXS:BBH:0612 0.6250  0.4997 -0.5001 0.000280 spinning_asym
# SXS:BBH:0614 0.5000  0.7500 -0.5000 0.000767 spinning_asym
# SXS:BBH:0615 0.5000  0.7500 -0.0001 0.000747 spinning_asym
# SXS:BBH:1142 0.8000 -0.0002 -0.0002 0.000082  non_spinning
# SXS:BBH:1143 0.8000 -0.0001 -0.0000 0.000088  non_spinning
# SXS:BBH:1147 0.4907  0.4301 -0.5100 0.000170 spinning_asym
# SXS:BBH:1166 0.5000  0.0000 -0.0000 0.000380  non_spinning
# SXS:BBH:1167 0.5000 -0.0001  0.0000 0.000362  non_spinning
# SXS:BBH:1178 0.3333 -0.0000  0.0000 0.000140  non_spinning
# SXS:BBH:1179 0.3333 -0.0000  0.0000 0.000037  non_spinning
# SXS:BBH:1221 0.3334 -0.0000 -0.0001 0.000005  non_spinning
# SXS:BBH:1222 0.5001  0.0000 -0.0001 0.000061  non_spinning
# SXS:BBH:1354 0.5458 -0.0002  0.0001 0.000053  non_spinning
# SXS:BBH:1387 0.3352  0.4699 -0.3599 0.000116 spinning_asym
# SXS:BBH:1412 0.6135  0.4000 -0.3000 0.000468 spinning_asym
# SXS:BBH:1414 0.5465 -0.5000  0.4000 0.000049 spinning_asym
# SXS:BBH:1446 0.3170 -0.7997  0.7765 0.000278 spinning_asym
# SXS:BBH:1458 0.2630 -0.0633  0.7994 0.000223 spinning_asym
# SXS:BBH:1461 0.3469 -0.4487 -0.7995 0.000641 spinning_asym
# SXS:BBH:1462 0.3801 -0.7999  0.5107 0.000662 spinning_asym
# SXS:BBH:1465 0.5855 -0.7876  0.7675 0.000128 spinning_asym
# SXS:BBH:1466 0.5272  0.6988 -0.7997 0.000230 spinning_asym
# SXS:BBH:1467 0.4484 -0.5624  0.7995 0.000370 spinning_asym
# SXS:BBH:1472 0.4227 -0.7998 -0.1160 0.000110 spinning_asym
# SXS:BBH:1478 0.5069  0.7998  0.1263 0.000676 spinning_asym
# SXS:BBH:1480 0.6459 -0.7998 -0.3096 0.000489 spinning_asym
# SXS:BBH:1483 0.3159  0.5607 -0.1941 0.000319 spinning_asym
# SXS:BBH:1484 0.3447 -0.5603  0.2974 0.000588 spinning_asym
# SXS:BBH:1485 0.3233  0.3527 -0.4007 0.000542 spinning_asym
# SXS:BBH:1486 0.2686  0.4272 -0.0336 0.000215 spinning_asym
# SXS:BBH:1488 0.6276 -0.3276  0.7472 0.000210 spinning_asym
# SXS:BBH:1489 0.2886  0.2998 -0.1717 0.000241 spinning_asym
# SXS:BBH:1491 0.6025  0.1993 -0.6965 0.000244 spinning_asym
# SXS:BBH:1906 0.2500  0.0001 -0.0001 0.000158  non_spinning
# SXS:BBH:1907 0.2501  0.0001  0.7994 0.000498 spinning_asym
# SXS:BBH:1911 0.2501 -0.0000 -0.7993 0.000040 spinning_asym
# SXS:BBH:1937 0.2500  0.4000 -0.0001 0.000393 spinning_asym
# SXS:BBH:1938 0.2501  0.3999  0.7995 0.000412 spinning_asym
# SXS:BBH:1942 0.2501  0.4000 -0.7994 0.000039 spinning_asym
# SXS:BBH:1961 0.2500 -0.4001  0.0000 0.000336 spinning_asym
# SXS:BBH:1962 0.2501 -0.4001  0.7996 0.000195 spinning_asym
# SXS:BBH:1966 0.2501 -0.4001 -0.7993 0.000109 spinning_asym
# SXS:BBH:2013 0.2500  0.0001  0.3999 0.000072 spinning_asym
# SXS:BBH:2036 0.2500  0.0000 -0.4000 0.000227 spinning_asym
# SXS:BBH:2110 0.5000 -0.5999 -0.0000 0.000062 spinning_asym
# SXS:BBH:2111 0.5000 -0.6000  0.5999 0.000041 spinning_asym
# SXS:BBH:2116 0.5000 -0.3000  0.3000 0.000039 spinning_asym
# SXS:BBH:2118 0.5000  0.0001 -0.6000 0.000106 spinning_asym
# SXS:BBH:2121 0.5000  0.0001  0.5998 0.000140 spinning_asym
# SXS:BBH:2123 0.5000  0.3001 -0.3000 0.000291 spinning_asym
# SXS:BBH:2128 0.5000  0.6001 -0.5998 0.000135 spinning_asym
# SXS:BBH:2129 0.5000  0.6000 -0.0001 0.000083 spinning_asym
# SXS:BBH:2134 0.3333 -0.6001 -0.0001 0.000027 spinning_asym
# SXS:BBH:2135 0.3334 -0.6000  0.5999 0.000132 spinning_asym
# SXS:BBH:2138 0.3333 -0.6001  0.3999 0.000151 spinning_asym
# SXS:BBH:2141 0.3334 -0.4001  0.5999 0.000084 spinning_asym
# SXS:BBH:2144 0.3333 -0.2999  0.3000 0.000139 spinning_asym
# SXS:BBH:2147 0.3333 -0.0000 -0.5999 0.000132 spinning_asym
# SXS:BBH:2150 0.3333  0.0001  0.6000 0.000098 spinning_asym
# SXS:BBH:2153 0.3333  0.3000 -0.3000 0.000106 spinning_asym
# SXS:BBH:2156 0.3334  0.4001 -0.5998 0.000288 spinning_asym
# SXS:BBH:2159 0.3333  0.6001 -0.5997 0.000086 spinning_asym
# SXS:BBH:2160 0.3333  0.6000 -0.3999 0.000204 spinning_asym
# SXS:BBH:2161 0.3333  0.5999  0.0000 0.000131 spinning_asym
# SXS:BBH:2265 0.3333  0.0000  0.0000 0.000072  non_spinning