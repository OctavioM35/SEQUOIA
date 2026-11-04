from typing import Any

from numpy import pi
def nnsur_convert(
    times,
    **kwargs) -> dict[str, Any]:
    out_pars = dict(
        q=kwargs['mass_ratio'],
        chiA0=kwargs['chi_1'],
        chiB0=kwargs['chi_2'],
        M=kwargs['total_mass'],
        dist_mpc=kwargs['luminosity_distance'],
        f_low=0,
        times=times,
        f_ref=None,
        freqs=None,
        inclination=pi-kwargs['theta_jn'],
        phi_ref=kwargs['phase'],
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