import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
from gwpy.timeseries import TimeSeries
from pycbc.detector import Detector
from scripts.surrogate.sur_utils import DANSur
from matplotlib.ticker import MultipleLocator


def synthetic_waveform(particular_event,output_synthetic):

    # ==========================================================
    # Parámetros sintéticos
    # ==========================================================

    synthetic_chirp_mass = 30.69       # masas solares
    synthetic_mass_ratio = 1/0.88      # m1/m2

    synthetic_gps_time = 1126259462.3047757
    # synthetic_geoscent_time = 1126259462.3047757
    synthetic_spin_1 = -0.05
    synthetic_spin_2 = -0.01

    synthetic_luminosity_distance = 473  # Mpc

    synthetic_psi = 1.45
    synthetic_declination = -1.19         # rad
    synthetic_ra = 2.03                   # rad
    synthetic_theta_jn = 2.70             # rad


    # ==========================================================
    # DANSur
    # ==========================================================

    surrogate = DANSur(device='cpu')

    q = synthetic_mass_ratio

    chi1 = [0, 0, synthetic_spin_1]
    chi2 = [0, 0, synthetic_spin_2]

    # Chirp mass en kg porque usamos units='mks'
    chirp_mass = synthetic_chirp_mass 

    # DANSur espera Mpc aquí
    dist_mpc = synthetic_luminosity_distance
    dt = 2 / 8192       # resolución temporal
    # ==========================================================
    # Tiempo
    # ==========================================================


    times = np.arange(-8, 8, dt)


    # ==========================================================
    # Generar waveform
    # ==========================================================

    waveform = surrogate(
        q,
        chi1,
        chi2,
        chirp_mass,
        dist_mpc,
        units="mks",
        f_low=0,
        times=np.arange(-8.0,8.0, dt),
        inclination=synthetic_theta_jn
    )

    # ==========================================================
    # Extraer waveform
    # ==========================================================
    times_out = np.squeeze(waveform[0])
    h = np.squeeze(waveform[1])
    print("waveform type:", type(waveform))
    print("waveform length:", len(waveform))

    print("times shape:", np.shape(waveform[0]))
    print("h shape:", np.shape(waveform[1]))

    print("times first/last:", waveform[0][0], waveform[0][-1])
    print("dt output:", waveform[0][1] - waveform[0][0])
    # numero de muestras necesarias: 65536
    # raise Exception


    # ==========================================================
    # Separar h+ y h×
    # ==========================================================

    # El surrogate devuelve:
    # h = h_plus - i*h_cross

    hp = h.real
    hc = -h.imag


    # ==========================================================
    # Tiempo relativo al merger
    # ==========================================================

    times_relative = times_out

    # ==========================================================
    # Tiempo GPS absoluto
    # ==========================================================

    # synthetic_gps_time corresponde al geocent_time,
    # es decir, al t = 0 del waveform.

    times_gps = synthetic_gps_time + times_out

    print("\nTiempo relativo:")
    print("  inicio:", times_relative[0])
    print("  merger:", 0.0)
    print("  final: ", times_relative[-1])

    print("\nTiempo GPS:")
    print("  inicio:", times_gps[0])
    print("  merger:", synthetic_gps_time)
    print("  final: ", times_gps[-1])


    # ==========================================================
    # Factores de antena de los detectores
    # ==========================================================

    detector_L1 = Detector("L1")
    detector_H1 = Detector("H1")

    fp_L1, fc_L1 = detector_L1.antenna_pattern(
        synthetic_ra,
        synthetic_declination,
        synthetic_psi,
        synthetic_gps_time
    )

    fp_H1, fc_H1 = detector_H1.antenna_pattern(
        synthetic_ra,
        synthetic_declination,
        synthetic_psi,
        synthetic_gps_time
    )

    print("\nAntenna patterns:")

    print("L1:")
    print("  F+ =", fp_L1)
    print("  Fx =", fc_L1)

    print("H1:")
    print("  F+ =", fp_H1)
    print("  Fx =", fc_H1)


    # ==========================================================
    # Strain visto por cada detector
    # ==========================================================

    strain_L1 = fp_L1 * hp + fc_L1 * hc
    strain_H1 = fp_H1 * hp + fc_H1 * hc


    print("\nStrain:")

    print(
        "L1:",
        np.min(strain_L1),
        np.max(strain_L1)
    )

    print(
        "H1:",
        np.min(strain_H1),
        np.max(strain_H1)
    )


    # ==========================================================
    # Plot L1
    # ==========================================================

    plt.figure(figsize=(10, 5))

    plt.plot(
        times_out,
        strain_L1,
        label="Synthetic strain L1"
    )

    plt.axvline(
        0,
        linestyle="--",
        label="Merger"
    )

    # Mostrar la zona alrededor de t = 7 s
    plt.xlim(6, 8)

    plt.xlabel("Time relative to merger [s]")
    plt.ylabel("Strain")
    plt.title("Synthetic gravitational-wave strain - L1")

    ax = plt.gca()
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    plt.savefig(
        "synthetic_strain_L1.png",
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()


    # ==========================================================
    # Guardar HDF5
    # ==========================================================

    event_dir = os.path.join(
        output_synthetic,
        particular_event
    )

    os.makedirs(event_dir, exist_ok=True)

    output_file = os.path.join(
        event_dir,
        particular_event + "_synthetic.hdf5"
    )


    with h5py.File(output_file, "w") as f:

        # ======================================================
        # Tiempos
        # ======================================================

        # Tiempo relativo al merger: -4 -> +4 s
        f.create_dataset(
            "times",
            data=np.asarray(times_relative, dtype=np.float64)
        )

        # Tiempo GPS absoluto
        f.create_dataset(
            "times_gps",
            data=np.asarray(times_gps, dtype=np.float64)
        )


        # ======================================================
        # Polarizaciones
        # ======================================================

        f.create_dataset(
            "h_plus",
            data=np.asarray(hp, dtype=np.float64)
        )

        f.create_dataset(
            "h_cross",
            data=np.asarray(hc, dtype=np.float64)
        )


        # ======================================================
        # Strain L1
        # ======================================================

        f.create_dataset(
            "strain_L1",
            data=np.asarray(strain_L1, dtype=np.float64)
        )


        # ======================================================
        # Strain H1
        # ======================================================

        f.create_dataset(
            "strain_H1",
            data=np.asarray(strain_H1, dtype=np.float64)
        )


        # ======================================================
        # Parámetros
        # ======================================================

        f.attrs["chirp_mass_msun"] = float(
            synthetic_chirp_mass
        )

        f.attrs["mass_ratio"] = float(
            synthetic_mass_ratio
        )

        f.attrs["luminosity_distance_mpc"] = float(
            synthetic_luminosity_distance
        )

        f.attrs["ra"] = float(
            synthetic_ra
        )

        f.attrs["dec"] = float(
            synthetic_declination
        )

        f.attrs["psi"] = float(
            synthetic_psi
        )

        f.attrs["theta_jn"] = float(
            synthetic_theta_jn
        )

        # Este es el GPS del merger/geocent time
        f.attrs["geocent_time"] = float(
            synthetic_gps_time
        )


        # ======================================================
        # Factores de antena
        # ======================================================

        f.attrs["Fplus_L1"] = float(
            fp_L1
        )

        f.attrs["Fcross_L1"] = float(
            fc_L1
        )

        f.attrs["Fplus_H1"] = float(
            fp_H1
        )

        f.attrs["Fcross_H1"] = float(
            fc_H1
        )


    # ==========================================================
    # Comprobar archivo
    # ==========================================================

    print("\nSynthetic waveform created successfully.")
    print(f"Saved to: {output_file}")

    with h5py.File(output_file, "r") as f:

        print("\nDatasets saved:")

        for key in f.keys():
            print(
                f"  {key}: "
                f"shape={f[key].shape}, "
                f"dtype={f[key].dtype}"
            )

        print("\nAttributes saved:")

        for key in f.attrs.keys():
            print(
                f"  {key}: {f.attrs[key]}"
            )