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
    # ---------------------------------------------------------
    # 1. Cargar waveform SXS
    # ---------------------------------------------------------
    sim = sxs.load(dicc["sxs"])
    h = sim.h


    # Dirección del observador
    theta = 0.0
    phi = 0.0
    ra = 0.0
    dec = 0.0
    psi = 0.0


    ht = h.evaluate(theta, phi)
    # LOS RESULTADOS SON INDEPENDIENTES DE LAS NORMALIZACIONES
    mass1 = (h.metadata["reference_mass1"])* 80   
    mass2=  (h.metadata["reference_mass2"])* 80   
    print('masas: ' , mass1, mass2)
    G = 6.67430e-11
    M_sun = 1.98847e30
    c = 299792458
    Mpc = 3.085677581e22

    total_mass = mass1 + mass2 
    luminosity_distance_mpc = 1000.0  

    # Expresión matemática dinámica
    factor = (2*G * (total_mass * M_sun)) / (c**2 * (luminosity_distance_mpc * Mpc))
    relaxation = h.metadata["relaxation_time"]


    # Recorte por tiempo, no por índice
    t_all = np.asarray(ht.t)
    mask = t_all >= relaxation


    t_all = np.asarray(ht.t)
    h_real = np.asarray(ht.real)
    h_imag = np.asarray(ht.imag)


    mask = t_all >= relaxation


    t_sxs = t_all[mask]


    h_plus = factor * h_real[mask]
    h_cross = -factor * h_imag[mask]
    total_mass = 40.0
    M_sun_time = 4.9254909476412675e-6
    time_conversion = total_mass * M_sun_time


    t_seconds = t_sxs * time_conversion


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


    geocent_time = 1000000000.0


    t_rel = t_uniform - t_uniform[-1]
    t_abs = geocent_time + t_rel


    dicc["start_time"] = t_abs[0]
    dicc["geocent_time"] = geocent_time





    # =========================================================
    # Respuesta de H1
    # =========================================================


    H1 = bilby.gw.detector.get_empty_interferometer("H1")


    # La respuesta de antena se evalúa en el tiempo geocéntrico
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


    print("F_plus  =", F_plus)
    print("F_cross =", F_cross)


    # =========================================================
    # Strain de H1
    # =========================================================


    h_H1 = (
        F_plus * h_plus_uniform
        + F_cross * h_cross_uniform
    )


    # =========================================================
    # TimeSeries con tiempo GPS absoluto
    # =========================================================


    ts = TimeSeries(
        h_H1,
        sample_rate=fs,
        t0=t_abs[0]
    )


    H1.set_strain_data_from_gwpy_timeseries(ts)


    print(f"len ts = {len(ts)}")
    print("H1 start time =", ts.t0.value)
    print("H1 end time   =", ts.times.value[-1])


    # =========================================================
    # Interferómetros
    # =========================================================


    ifos = InterferometerList([H1])


    # Guardamos también la información de la inyección
    dicc["signal"] = [
        ts,
        h_cross_uniform,
        h_plus_uniform
    ]


    dicc["geocent_time"] = geocent_time
    dicc["ra"] = ra
    dicc["dec"] = dec
    dicc["psi"] = psi


    return ifos 