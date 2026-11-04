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
    return out_pars