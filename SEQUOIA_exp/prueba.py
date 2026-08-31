
import matplotlib.pyplot as plt
import numpy as np
import gwsurrogate as gws
from astropy.cosmology import Planck18, z_at_value
import time
import astropy.units as u
from scripts.surrogate.sur_utils import DANSur


def masa_total(masa_chirp_det, q):
    q_nrh = 1.0/q

    # -----------------------------------------
    # Redshift a partir de la luminosidad
    # -----------------------------------------
    z = z_at_value(
        Planck18.luminosity_distance,
        dist * u.Mpc
    )

    z = float(z)

    # -----------------------------------------
    # Chirp mass: detector -> source frame
    # -----------------------------------------
    masa_chirp_source = masa_chirp_det / (1 + z)

    # -----------------------------------------
    # Masa total: source frame
    # -----------------------------------------
    M = (
        masa_chirp_source
        * (1 + q_nrh)**(6 / 5)
        / q_nrh**(3 / 5)
    )
    return M , q_nrh





#CONFIGURACIÓN DE PARÁMETROS
chiA = [0, 0, 0]
chiB = [0, 0, 0]

dt = 1 / 4096
f_low = 20
f_ref = 20
dist = 4457  # Mpc

q = 0.6
masa_chirp_det = 54  # Msun, detector frame



M , q_nrh = masa_total(masa_chirp_det,q)
# print("================================")
# print('PARAMETROS CARGADOS')
# print("dist =", dist, "Mpc")
# print("q original =", q)
# print("Mchirp detector =", masa_chirp_det, "Msun")
# print("Mtotal fuente =", M, "Msun")
# print("================================")

# print('LLAMADA CON NRH')
t0 = time.perf_counter()
nrh=gws.LoadSurrogate("NRSur7dq4")
import numpy as np

# ============================================================
# Parámetros
# ============================================================

q = 1.3837116999873722

chiA0 = np.array([
    0.31219747572788453,
    0.0,
    0.011489638382098473
])

chiB0 = np.array([
    -0.014640975258558234,
    -0.41907125223060143,
    -0.1586028982064239
])

M_solar = 47.55362873081309
dist_mpc = 9945.467991634205

f_low_hz = 20.0
f_ref_hz = 20.0

dt_seconds = 0.00048828125

# ============================================================
# Conversión a unidades adimensionales
# ============================================================

# 1 M_sun en segundos (G=c=1)
M_sun_seconds = 4.92549095e-6

# Masa total en segundos
M_seconds = M_solar * M_sun_seconds

# Frecuencias adimensionales M*f
f_low = M_seconds * f_low_hz
f_ref = M_seconds * f_ref_hz

# Paso temporal adimensional: dt/M
dt = dt_seconds / M_seconds

# ============================================================
# Mostrar parámetros
# ============================================================

print("========== PARÁMETROS NRSur DIMENSIONLESS ==========")
print("q       =", q)
print("chiA0   =", chiA0)
print("chiB0   =", chiB0)
print("M       =", M_solar, "Msun")
print("M       =", M_seconds, "s")
print("dist    =", dist_mpc, "Mpc")
print("f_low   =", f_low, "(M*f)")
print("f_ref   =", f_ref, "(M*f)")
print("dt      =", dt, "(dt/M)")
print("=====================================================")

# ============================================================
# NRSur7dq4
# ============================================================

domain, h, dynamics = nrh(
    q=q,
    chiA0=chiA0,
    chiB0=chiB0,
    f_low=f_low,
    f_ref=f_ref,
    dt=dt,
    units='dimensionless'
)

print("domain shape:", np.shape(domain))
print("h shape:", np.shape(h))
print("dynamics type:", type(dynamics))



print(f"Tiempo requerido por NRH: {time.perf_counter() - t0:.6f} s")

print(len(domain))
print(len(h))
print(h)
plt.figure()
plt.plot(np.imag(h[(2, 2)]))
plt.savefig('prueba nrsur.png')


# print('LLAMADA CON DANSUR')
# surrogate = DANSur(device='cpu')

# t1 = time.perf_counter()
# waveform = surrogate(q=q_nrh, chiA0=chiA, chiB0=chiB)
# print(f"Tiempo requerido por DANSur: {time.perf_counter() - t0:.6f} s")

# plt.figure()
# plt.plot(waveform)
# plt.savefig('DANSur.png')


# q = 1.41
# chiA = [0, 0, 0.5]
# chiB = [0, 0, -0.7]
# dt =1/4096 
# f_low=20
# f_ref=20     
# dist=13695
# M=166