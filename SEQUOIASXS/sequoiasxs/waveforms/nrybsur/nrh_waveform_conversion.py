from math import pi
def nnsur_convert(
    # dt,
    times,
    mass_ratio=None,
    total_mass=None,
    chi_1=0.1,
    chi_2=0.1,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    luminosity_distance=300.0,
    theta_jn=0.0,
    psi=0,
    phase=0,
    geocent_time=1126259642.413,
    ra=0,
    dec=0, 
    **kwargs):
    r"""
    This example only creates a linearly polarised signal with only plus
    polarisation.
    .. math::

        h_{\plus}(t) =
            \Theta(t - t_{0}) A
            e^{-(t - t_{0}) / \tau}
            \sin \left( 2 \pi f t + \phi \right)

    Parameters
    ----------
    mass_1=40.0,
    mass_2=20.0,
    chi_1=0.1,
    chi_2=0.1,
    tilt_1=0.0,
    tilt_2=0.0,
    phi_12=0.0,
    phi_jl=0.0,
    luminosity_distance=300.0,
    theta_jn=0.0,
    psi=0,
    phase=0,
    geocent_time=1126259642.413,
    ra=0,
    dec=0,
)

    Returns
    -------
    dict:
        A dictionary containing "plus" and "cross" entries.

    """
    # print(time)
    # chi_1 = chi_1
    # chi_2 = chi_2
    q_val = 1.0 / mass_ratio if mass_ratio <= 1.0 else float(mass_ratio)
    if total_mass is None and kwargs.get('chirp_mass') is not None:
        q_norm = mass_ratio if mass_ratio <= 1.0 else 1.0 / mass_ratio
        total_mass = float(kwargs['chirp_mass']) * ((1.0 + q_norm) ** 1.2) / (q_norm ** 0.6)
    out_pars = dict(
        q=q_val, 
        chiA0=chi_1,
        chiB0=chi_2,
        M=total_mass,
        dist_mpc=luminosity_distance,
        f_low=0,
        times=times,
        f_ref=None,
        freqs=None,
        inclination=float(theta_jn),
        phi_ref=float(phase),
        units='mks',
    )
    return out_pars