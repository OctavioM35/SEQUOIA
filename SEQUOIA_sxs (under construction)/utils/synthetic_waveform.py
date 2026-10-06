import numpy as np
import bilby
import sxs
from config import sxs_simulation
from bilby.gw.detector import InterferometerList
from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from config import dicc
from bilby.gw.detector.psd import PowerSpectralDensity
from parameter_estimation_pipeline.mismatch import *


def synthetic_waveform(dicc):
    """
        Generate synthetic gravitational-wave data from an SXS numerical-relativity
    waveform.

        The function loads an SXS simulation, evaluates the waveform for a chosen
    observer direction, rescales it to a physical total mass and luminosity
    distance, projects the plus and cross polarizations onto the H1 and L1
    detectors, adds simulated detector noise, and stores the corresponding
    injection parameters in the configuration dictionary.

         The function returns an InterferometerList containing the H1 and L1
    interferometers with the noisy injected data.

    INPUTS:
    -------
    
    -dicc: The created dicctionay in config.py 

    """

    sim = sxs.load(dicc["sxs"])
    h = sim.h


    # Observer direction and source position
    theta = 1
    phi = 0.0
    ra = 0.0
    dec = 0.0
    psi = 0.0


    ht = h.evaluate(theta, phi)


    
    # Scailing normalizations

    luminosity_distance_mpc = 1000
    mass_scailing = 100



    #Conversion units
    G = 6.67430e-11
    M_sun = 1.98847e30
    c = 299792458
    Mpc = 3.085677581e22

    mass1 = (h.metadata["reference_mass1"])* mass_scailing 
    mass2=  (h.metadata["reference_mass2"])* mass_scailing
    total_mass = mass1 + mass2 



    # Amplitude scaling factor
    factor = (2*G * (total_mass * M_sun)) / (c**2 * (luminosity_distance_mpc * Mpc))
    relaxation = h.metadata["relaxation_time"]

    # Select the waveform data after the relaxation time

    t_all = np.asarray(ht.t)
    mask = t_all >= relaxation
    t_all = np.asarray(ht.t)
    h_real = np.asarray(ht.real)
    h_imag = np.asarray(ht.imag)

    mask = t_all >= relaxation
    t_sxs = t_all[mask]

    h_plus = factor * h_real[mask]
    h_cross = -factor * h_imag[mask]

    # Convert the SXS time coordinate to seconds
    M_sun_time = 4.9254909476412675e-6
    time_conversion = total_mass * M_sun_time
    t_seconds = t_sxs * time_conversion


    # Interpolate the waveform onto a uniform sampling grid
    fs = 4096


    t_uniform = np.arange(
        t_seconds[0],
        t_seconds[-1],
        1.0 / fs
    )


    h_plus_uniform = np.interp(
        t_uniform,
        t_seconds,
        h_plus
    )


    h_cross_uniform = np.interp(
        t_uniform,
        t_seconds,
        h_cross
    )

    # Define the absolute GPS time
    geocent_time = 1000000000.0
    dicc["geocent_time"] = geocent_time


    t_rel = t_uniform - t_uniform[-1]
    t_abs = geocent_time + t_rel


    dicc["start_time"] = t_abs[0]
    dicc["geocent_time"] = geocent_time


    # Initialize the H1 and L1 interferometers
    H1 = bilby.gw.detector.get_empty_interferometer("H1")
    psd = bilby.gw.detector.PowerSpectralDensity(frequency_array = np.linspace(15, 4096, 1000), psd_array = np.ones(1000)*1e-46)
    H1.power_spectral_density=psd
    
    L1 = bilby.gw.detector.get_empty_interferometer("L1")
    psd = bilby.gw.detector.PowerSpectralDensity(frequency_array = np.linspace(15, 4096, 1000), psd_array = np.ones(1000)*1e-46)
    L1.power_spectral_density=psd


    # Evaluate the detector antenna responses
    F_plus = H1.antenna_response(
        ra=ra,
        dec=dec,
        time=geocent_time,
        psi=psi,
        mode="plus"
    )


    F_cross = H1.antenna_response(
        ra=ra,
        dec=dec,
        time=geocent_time,
        psi=psi,
        mode="cross"
    )

    F_plus_l1 = L1.antenna_response(
        ra=ra,
        dec=dec,
        time=geocent_time,
        psi=psi,
        mode="plus"
    )


    F_cross_l1 = L1.antenna_response(
        ra=ra,
        dec=dec,
        time=geocent_time,
        psi=psi,
        mode="cross"
    )

    # Project the polarizations onto the detectors
    h_H1 = (
        F_plus * h_plus_uniform
        + F_cross * h_cross_uniform
    )

    h_L1 = (
        F_plus_l1 * h_plus_uniform
        + F_cross_l1 * h_cross_uniform
    )


    ts = TimeSeries(
        h_H1,
        sample_rate=fs,
        t0=t_abs[0]
    )

    ts_l1 = TimeSeries(
        h_L1,
        sample_rate=fs,
        t0=t_abs[0]
    )


    # Generate and add detector noise for H1
    H1.set_strain_data_from_power_spectral_density(fs, ts.duration.value)
    noise=H1.time_domain_strain
    inj=noise+ts
    H1.set_strain_data_from_gwpy_timeseries(inj)


    # Generate and add detector noise for L1
    L1.set_strain_data_from_power_spectral_density(fs, ts_l1.duration.value)
    noise=L1.time_domain_strain
    inj=noise+ts
    L1.set_strain_data_from_gwpy_timeseries(inj)

    print(f"len ts = {len(ts)}")
    print("H1 start time =", ts.t0.value)
    print("H1 end time   =", ts.times.value[-1])

    ifos = InterferometerList([H1,L1])

    #Store inyection values and signal
    dicc["signal"] = [
        ts,
        h_cross_uniform,
        h_plus_uniform
    ]
    chirp_mass = (
        (mass1 * mass2)**(3.0 / 5.0)
        / (mass1 + mass2)**(1.0 / 5.0)
    )

    chi_1 = h.metadata["reference_dimensionless_spin1"][2]
    chi_2 = h.metadata["reference_dimensionless_spin2"][2]
    inyection_parameters = {'mass_ratio': mass2/mass1,
                            'chirp_mass' : chirp_mass,
                            'luminosity_distance': luminosity_distance_mpc,
                            "chi_1" : chi_1,
                            "chi_2" : chi_2,
                            }

    dicc["inyection"] = inyection_parameters

    return ifos 