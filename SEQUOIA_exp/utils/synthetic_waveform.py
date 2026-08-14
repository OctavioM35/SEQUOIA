import sys
import os
import numpy as np
import bilby
import h5py

from nessai.utils.logging import configure_logger
import numpy as np

from scripts.surrogate.sur_utils import DANSur
from waveform_conversion import nnsur_convert
from gwpy.frequencyseries import FrequencySeries
from bilby.gw.prior import *

from bilby.gw.waveform_generator import WaveformGenerator
from bilby.gw.likelihood import GravitationalWaveTransient
from bilby.gw.detector import InterferometerList

from waveform_generator import waveform_generator
from matplotlib import pyplot as plt
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine, Cosine
import ast

# ============================================================
# CONFIGURACIÓN DEL LOGGER
# ============================================================
#
# Mostrar WARNING y ERROR.
#
# Esto elimina mensajes INFO como:
#
# [LIKELIHOOD] Evaluation #14287 | ...
#
# pero mantiene los WARNING y ERROR.
# ============================================================

configure_logger(log_level="WARNING")


def synthetic_waveform(dicc,output_synthetic,folder,zenodo_file):


    # ==========================================================
    # CONFIGURACIÓN
    # ==========================================================

    output_synthetic = (
        "/home/octaio-m/PycharmProjects/"
        "PythonProject/TFG/DANSur_22-master/"
        "sintetico"
    )

    duration = float(dicc["duration"])

    sampling_frequency = float(
        dicc["sampling-frequency"]
    )

    # ==========================================================
    # CREAR DIRECTORIO DE SALIDA
    # ==========================================================
    outdir = os.path.join(output_synthetic, f"synthetic_outdir_NN_{folder}")
    os.makedirs(
        outdir,
        exist_ok=True
    )

    # ==========================================================
    # INJECTION PARAMETERS
    # ==========================================================
    with h5py.File(zenodo_file, "r") as f:
        key0 = list(f.keys())[0]
        injection_parameters = {

                "chirp_mass": np.median(f[key0]["posterior_samples"]["chirp_mass"]),

                "mass_ratio": np.median(f[key0]["posterior_samples"]["mass_ratio"]),

                "chi_1": 0.0,
                "chi_2":0.0,

                "tilt_1": 0.0,
                "tilt_2": 0.0,

                "phi_12": 0.0,
                "phi_jl": 0.0,

                "luminosity_distance": np.median(f[key0]["posterior_samples"]["luminosity_distance"]),

                "theta_jn":np.median(f[key0]["posterior_samples"]["theta_jn"]),

                "psi": np.median(f[key0]["posterior_samples"]["psi"]),

                "phase": 0.0,

                "geocent_time": np.median(f[key0]["posterior_samples"]["geocent_time"]),

                "ra":np.median(f[key0]["posterior_samples"]["ra"]),
                "dec": np.median(f[key0]["posterior_samples"]["dec"]),
            }

    # ==========================================================
    # PRIORS
    # ==========================================================

    priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)
    with h5py.File(zenodo_file, "r") as f:
            key0 = list(f.keys())[0]
            try:
                    geocent_time_prior = f[key0]['priors']['analytic']['geocent_time'][:][0]
                    print('Prior: geocent')

            except:
                    try:
                        geocent_time_prior = f[key0]['priors']['analytic']['L1_time'][:][0]
                    except:
                        try:
                            geocent_time_prior = f[key0]['priors']['analytic']['H1_time'][:][0]
                        except:
                            return None, None
            priors['geocent_time'] = eval(geocent_time_prior)



            priors["luminosity_distance"] = eval(f[key0]["priors"]["analytic"]["luminosity_distance"][:][0] )

            priors["theta_jn"] = eval(f[key0]["priors"]["analytic"]["theta_jn"][:][0] )
            priors["psi"] = eval(f[key0]["priors"]["analytic"]["psi"][:][0] )


            priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
            try:
                priors["dec"] = eval(f[key0]["priors"]["analytic"]["dec"][:][0])
                priors["ra"] = eval(f[key0]["priors"]["analytic"]["ra"][:][0])
                priors["chirp_mass"] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )


                priors["mass_ratio"] = eval(f[key0]["priors"]["analytic"]["mass_ratio"][:][0] )
            except KeyError:
                priors["azimuth"] = eval(f[key0]["priors"]["analytic"]["azimuth"][:][0])
                priors["zenith"] = eval(f[key0]["priors"]["analytic"]["zenith"][:][0])
            priors['chi_1'] = bilby.gw.prior.AlignedSpin(
                                    a_prior=bilby.core.prior.Uniform(minimum=-0.8, maximum=0.8, name=None, 
                                                                    latex_label=None, unit=None, boundary=None), 
                                    z_prior=bilby.core.prior.Uniform(minimum=-1, maximum=1, name=None, 
                                                                    latex_label=None, unit=None, boundary=None), 
                                    name='chi_1', latex_label='$\\chi_1$', unit=None, boundary=None, 
                                    minimum=-0.8, maximum=0.8)
            priors['chi_2'] = bilby.gw.prior.AlignedSpin(
                                    a_prior=bilby.core.prior.Uniform(minimum=-0.8, maximum=0.8, name=None, 
                                                                    latex_label=None, unit=None, boundary=None), 
                                    z_prior=bilby.core.prior.Uniform(minimum=-1, maximum=1, name=None, 
                                                                    latex_label=None, unit=None, boundary=None), 
                                    name='chi_2', latex_label='$\\chi_2$', unit=None, boundary=None, 
                                    minimum=-0.8, maximum=0.8)

    # ==========================================================
    # TARG_KEYS
    # ==========================================================

    targ_keys = [
        "chirp_mass",
        "mass_ratio",
        "chi_1",
        "chi_2",
        "tilt_1",
        "tilt_2",
        "phi_12",
        "phi_jl",
        "luminosity_distance",
        "theta_jn",
        "psi",
        "phase",
        "geocent_time",
        "ra",
        "dec",
    ]

    print()
    print("=" * 70)
    print("TARG KEYS")
    print("=" * 70)

    print(targ_keys)

    # ==========================================================
    # COMPROBAR INJECTION PARAMETERS
    # ==========================================================

    missing = [
        key for key in targ_keys
        if key not in injection_parameters
    ]

    if missing:
        raise RuntimeError(
            f"Faltan parámetros en injection_parameters: {missing}"
        )

    # ==========================================================
    # COMPROBAR PRIORS
    # ==========================================================

    missing_priors = [
        key for key in targ_keys
        if key not in priors
    ]

    if missing_priors:

        print()
        print("WARNING:")

        print(
            "Los siguientes targ_keys no tienen prior explícito:"
        )

        print(missing_priors)

    # ==========================================================
    # WAVEFORM GENERATOR
    # ==========================================================

    generated_waveform = waveform_generator(
        dicc,
        targ_keys
    )

    print()
    print("=" * 70)
    print("WAVEFORM GENERATOR")
    print("=" * 70)

    print(generated_waveform)

    print()
    print("Duration =", duration)

    print(
        "Sampling frequency =",
        sampling_frequency
    )

    print(
        "Start time =",
        generated_waveform.start_time
    )

    print(
        "Geocent time =",
        injection_parameters["geocent_time"]
    )

    print(
        "Time from start to geocent =",
        injection_parameters["geocent_time"]
        - generated_waveform.start_time
    )

    # ==========================================================
    # INTERFERÓMETRO
    # ==========================================================

    ifos = InterferometerList(["L1"])

    # ==========================================================
    # DATOS DEL DETECTOR
    # ==========================================================


    with h5py.File(zenodo_file) as file:
        ap_key = list(file.keys())[0]
        fs, data = file[f'{ap_key}/psds/L1'][()].T
        ifos_asd = FrequencySeries(data, frequencies=fs)
        ifos.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(ifos_asd.frequencies, ifos_asd.value)

    ifos.set_strain_data_from_power_spectral_densities(

        duration=duration,

        sampling_frequency=sampling_frequency,

        start_time=generated_waveform.start_time
    )
    # ==========================================================
    # INYECTAR SEÑAL
    # ==========================================================

    print()
    print("=" * 70)
    print("INJECTING SIGNAL")
    print("=" * 70)

    ifos.inject_signal(
        waveform_generator=generated_waveform,
        parameters=injection_parameters
    )

    print("Injection completed.")

    # ==========================================================
    # RESPUESTA DEL DETECTOR
    # ==========================================================
    polarizations = (
        generated_waveform.frequency_domain_strain(
            injection_parameters
        )
    )

    for ifo in ifos:

        detector_response = (
            ifo.get_detector_response(
                polarizations,
                injection_parameters
            )
        )

        print()
        print("=" * 70)
        print(
            f"DETECTOR RESPONSE: {ifo.name}"
        )
        print("=" * 70)

        print(
            "shape:",
            detector_response.shape
        )

        print(
            "max |strain|:",
            np.max(
                np.abs(detector_response)
            )
        )

    # ==========================================================
    # LIKELIHOOD
    # ==========================================================

    likelihood = GravitationalWaveTransient(
        interferometers=ifos,
        waveform_generator=generated_waveform,
        priors=priors,
        distance_marginalization=False,
        phase_marginalization=dicc["phase-marginalization"],
        time_marginalization=dicc["time-marginalization"],
        jitter_time=dicc["jitter-time"]
    )

    # # ==========================================================
    # # TRUE PARAMETERS
    # # ==========================================================

    # likelihood.parameters.update(
    #     injection_parameters
    # )

    # # ==========================================================
    # # TRUE LIKELIHOOD
    # # ==========================================================

    # logL_true = likelihood.log_likelihood()

    # # ==========================================================
    # # TEST CHIRP MASS
    # # ==========================================================

    # test_parameters = (
    #     injection_parameters.copy()
    # )

    # test_parameters["chirp_mass"] += 1.0

    # likelihood.parameters.update(
    #     test_parameters
    # )

    # logL_mc = likelihood.log_likelihood()

    # print()
    # print("=" * 70)
    # print("TEST CHIRP MASS")
    # print("=" * 70)

    # print("logL(TRUE) =", logL_true)

    # print("logL(Mc+1) =", logL_mc)

    # print(
    #     "Delta logL =",
    #     logL_true - logL_mc
    # )

    # # ==========================================================
    # # TEST DISTANCE
    # # ==========================================================

    # likelihood.parameters.update(
    #     injection_parameters
    # )

    # test_parameters = (
    #     injection_parameters.copy()
    # )

    # test_parameters[
    #     "luminosity_distance"
    # ] *= 2

    # likelihood.parameters.update(
    #     test_parameters
    # )

    # logL_distance = (
    #     likelihood.log_likelihood()
    # )

    # print()
    # print("=" * 70)
    # print("TEST DISTANCE")
    # print("=" * 70)

    # print("logL(TRUE) =", logL_true)

    # print(
    #     "logL(distance x2) =",
    #     logL_distance
    # )

    # print(
    #     "Delta logL =",
    #     logL_true - logL_distance
    # )

    # # ==========================================================
    # # TEST TIME
    # # ==========================================================

    # likelihood.parameters.update(
    #     injection_parameters
    # )

    # test_parameters = (
    #     injection_parameters.copy()
    # )

    # test_parameters[
    #     "geocent_time"
    # ] += 0.1

    # likelihood.parameters.update(
    #     test_parameters
    # )

    # logL_time = (
    #     likelihood.log_likelihood()
    # )

    # # ==========================================================
    # # TEST PHASE
    # # ==========================================================

    # likelihood.parameters.update(
    #     injection_parameters
    # )

    # test_parameters = (
    #     injection_parameters.copy()
    # )

    # test_parameters["phase"] += (
    #     np.pi / 2
    # )

    # likelihood.parameters.update(
    #     test_parameters
    # )

    # logL_phase = (
    #     likelihood.log_likelihood()
    # )

    # print()
    # print("=" * 70)
    # print("TEST PHASE")
    # print("=" * 70)

    # print("logL(TRUE) =", logL_true)

    # print(
    #     "logL(phase+pi/2) =",
    #     logL_phase
    # )

    # print(
    #     "Delta logL =",
    #     logL_true - logL_phase
    # )

    # ==========================================================
    # RESTAURAR TRUE PARAMETERS
    # ==========================================================

    likelihood.parameters.update(
        injection_parameters
    )

    # ==========================================================
    # RESUMEN
    # ==========================================================

    # print()
    # print("=" * 70)
    # print("SUMMARY")
    # print("=" * 70)

    # print(
    #     "logL TRUE         =",
    #     logL_true
    # )

    # print(
    #     "logL time + 0.1   =",
    #     logL_time
    # )

    # print(
    #     "logL Mc + 1       =",
    #     logL_mc
    # )

    # print(
    #     "logL distance x2  =",
    #     logL_distance
    # )

    # print(
    #     "logL phase + pi/2 =",
    #     logL_phase


    # ==========================================================
    # NESSAI
    # ==========================================================

    print()
    print("=" * 70)
    print("STARTING NESSAI")

    result = bilby.run_sampler(
        likelihood=likelihood,
        priors=priors,
        sampler="nessai",
        use_ratio=False,
        flow_class="gwflowproposal",
        npoints=200,
        naccept=10,
        injection_parameters=(
            injection_parameters),
        resume=True,
        outdir=outdir,
        label="DANSur",
        npool=1,
        stopping=1
    )

    # ==========================================================
    # RESULTADOS
    # ==========================================================

    print()
    print("=" * 70)
    print("SAMPLING FINISHED")
    print("=" * 70)

    print(result.posterior)

    # ==========================================================
    # ARCHIVOS GENERADOS
    # ==========================================================

    print()
    print("=" * 70)
    print("FILES GENERATED")
    print("=" * 70)

    for filename in sorted(
        os.listdir(output_synthetic)
    ):
        print(filename)

    # ==========================================================
    # CORNER
    # ==========================================================
    plt.figure()
    maxll_params = dict(result.posterior[result.posterior.log_likelihood == result.posterior.log_likelihood.max()].iloc[0])
    plt.plot(generated_waveform.time_array, generated_waveform.time_domain_strain(maxll_params)['plus'], '--')

    plt.title("Reconstructed waveform", fontsize=14)
    plt.xlabel("Time [s]", fontsize=12)
    plt.ylabel("h strain", fontsize=12)
    plt.savefig(os.path.join(outdir, 'reconstructed_waveform.png' ))