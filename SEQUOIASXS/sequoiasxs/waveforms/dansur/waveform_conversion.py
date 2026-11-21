from typing import Any
import numpy as np


def nnsur_convert(
    times,
    **kwargs) -> dict[str, Any]:
    theta_jn = float(kwargs.get('theta_jn', 0.0))
    phi_ref = float(kwargs.get('phase', kwargs.get('phi_ref', 0.0)))
    q_val = float(kwargs.get('mass_ratio', kwargs.get('mass_2', 1.0) / kwargs.get('mass_1', 1.0)))
    if q_val > 1.0:
        q_val = 1.0 / q_val

    if 'total_mass' in kwargs:
        M = float(kwargs['total_mass'])
    elif 'chirp_mass' in kwargs:
        mc = float(kwargs['chirp_mass'])
        M = mc * ((1.0 + q_val) ** 1.2) / (q_val ** 0.6)
    else:
        M = float(kwargs.get('mass_1', 35.0)) + float(kwargs.get('mass_2', 35.0))

    out_pars = dict(
        q=q_val,
        chiA0=float(kwargs.get('chi_1', 0.0)),
        chiB0=float(kwargs.get('chi_2', 0.0)),
        M=M,
        dist_mpc=float(kwargs.get('luminosity_distance', 100.0)),
        f_low=0,
        times=times,
        f_ref=None,
        freqs=None,
        inclination=theta_jn,
        phi_ref=phi_ref,
        units='mks',
    )
    """
    Convert the sampled Bilby parameters into the keyword arguments
    expected by the surrogate waveform evaluation (gwsurrogate-style call).

    The function renames and rearranges the parameters of the inference
    (mass ratio, aligned spins, total mass, distance, inclination and
    phase) into the dictionary required by the surrogate. The reference
    frequency and the frequency array are set to None and the low
    frequency to 0, so the waveform is evaluated on the given time array
    and returned in SI (mks) units.

    INPUTS:
    -------

    -times: Time array (in the units expected by the surrogate) at which
    the waveform is evaluated.

    -**kwargs: Sampled parameters. The following keys are required:
        - mass_ratio: Mass ratio q of the binary.
        - chi_1, chi_2: Dimensionless spin components of the two black
          holes along the orbital angular momentum.
        - total_mass: Total mass of the binary (solar masses).
        - luminosity_distance: Luminosity distance (Mpc).
        - theta_jn: Angle between the total angular momentum and the
          line of sight (rad).
        - phase: Reference orbital phase (rad).

    OUTPUTS:
    --------

    -Dictionary with the keys q, chiA0, chiB0, M, dist_mpc, f_low, times,
    f_ref, freqs, inclination, phi_ref and units, ready to be passed to
    the surrogate evaluation function. The inclination is computed as
    pi - theta_jn.
    """
    return out_pars