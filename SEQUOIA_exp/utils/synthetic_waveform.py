import sys
import os
import numpy as np
import bilby
import h5py

from nessai.utils.logging import configure_logger
import numpy as np

from scripts.surrogate.sur_utils import DANSur
from gwpy.frequencyseries import FrequencySeries
from bilby.gw.prior import *

from bilby.gw.waveform_generator import WaveformGenerator
from bilby.gw.likelihood import GravitationalWaveTransient
from bilby.gw.detector import InterferometerList

from matplotlib import pyplot as plt
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine, Cosine
import ast


def synthetic_waveform(dicc,zenodo_file):

    # ==========================================================
    # INJECTION PARAMETERS
    # ==========================================================
    if dicc["automatic"] == True:
        with h5py.File(zenodo_file, "r") as f:
            key0 = list(f.keys())[0]
            injection_parameters = {

                    "chirp_mass": np.median(f[key0]["posterior_samples"]["chirp_mass"]),

                    "mass_ratio": np.median(f[key0]["posterior_samples"]["mass_ratio"]),

                    "chi_1": np.median(f[key0]["posterior_samples"]["a_1"]),
                    "chi_2":np.median(f[key0]["posterior_samples"]["a_2"]),

                    "tilt_1": 0.0,
                    "tilt_2": 0.0,

                    "phi_12": 0.0,
                    "phi_jl": 0.0,

                    "luminosity_distance": np.median(f[key0]["posterior_samples"]["luminosity_distance"]),

                    "theta_jn":np.median(f[key0]["posterior_samples"]["theta_jn"]),

                    "psi": np.median(f[key0]["posterior_samples"]["psi"]),

                    "phase": np.median(f[key0]["posterior_samples"]["phase"]),

                    "geocent_time": np.median(f[key0]["posterior_samples"]["geocent_time"]),

                    "ra":np.median(f[key0]["posterior_samples"]["ra"]),
                    "dec": np.median(f[key0]["posterior_samples"]["dec"]),
                }

    
            dicc["injection_parameters"] = injection_parameters

    # ==========================================================
    # INTERFEROMETER
    # ==========================================================

    ifos = InterferometerList(["L1"])

    with h5py.File(zenodo_file) as file:
        ap_key = list(file.keys())[0]
        fs, data = file[f'{ap_key}/psds/L1'][()].T
        ifos_asd = FrequencySeries(data, frequencies=fs)
        ifos.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(ifos_asd.frequencies, ifos_asd.value)

    ifos.set_strain_data_from_power_spectral_densities(

        duration=dicc["duration"],

        sampling_frequency=dicc["sampling-frequency"],

        start_time=dicc["injection_parameters"]["geocent_time"]
    )
    return ifos