from math import pi

import numpy as np
def nrsur_convert(
    # dt,
    times,
    mass_ratio=None,
    total_mass=None,
    a_1=0.1,
    a_2=0.1,
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
    chiA0 = [ a_1 * np.sin(tilt_1), 0.0, a_1 * np.cos(tilt_1) ]
    chiB0 = [ a_2 * np.sin(tilt_2) * np.cos(phi_12), a_2 * np.sin(tilt_2) * np.sin(phi_12), a_2 * np.cos(tilt_2) ]   
    # chi_1 = chi_1
    # chi_2 = chi_2
    a_1=a_1,
    a_2=a_2,
    out_pars = dict(
    q = 1/mass_ratio, 
    chiA0 = chiA0,
    chiB0 = chiB0,
    M = total_mass,
    dist_mpc=luminosity_distance,
    f_low=0,
    # dt=dt, 
    times = times,
    f_ref=None,
    freqs=None,
    inclination= pi - theta_jn,
    phi_ref = phase,
    units='mks',
    )
    return out_pars