import bilby
import numpy as np


def load_custom_priors(dicc):
    """
    Establish custom priors and targ keys for the PE

    INPUTS:
    ------
    -dicc: The created dictionary in config.py 

    Returns
    ------
    -targ_keys: Names of the source parameters used by the waveform generator.
    -priors: Prior distributions for the parameters sampled during the Bayesian
             inference.
    """
    print('Loading LVK priors...')

    mc_inj = float(dicc.get("chirp_mass", dicc.get("inyection", {}).get("chirp_mass", 30.0)))
    tc_inj = float(dicc.get("geocent_time", 1126259462.0))
    dist_inj = float(dicc.get("inyection", {}).get("luminosity_distance", 1000.0))

    if dicc.get("precession_model", False) is False:
        # Non-precessing (aligned spin) model
        priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)

        priors['chirp_mass'] = bilby.gw.prior.UniformInComponentsChirpMass(
            minimum=10,
            maximum=100,
            name='chirp_mass',
            latex_label=r'$\mathcal{M}$'
        )
        priors['mass_ratio'] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=0.125, maximum=1.0, name='mass_ratio')
        # priors['chi_1'] = bilby.gw.prior.AlignedSpin(
        #                                     a_prior=bilby.core.prior.Uniform(minimum=0, maximum=0.8, name=None, 
        #                                                                     latex_label=None, unit=None, boundary=None), 
        #                                     z_prior=bilby.core.prior.Uniform(minimum=-1, maximum=1, name=None, 
        #                                                                     latex_label=None, unit=None, boundary=None), 
        #                                     name='chi_1', latex_label='$\\chi_1$', unit=None, boundary=None, 
        #                                     minimum=-0.8, maximum=0.8)
        # priors['chi_2'] = bilby.gw.prior.AlignedSpin(
        #                                     a_prior=bilby.core.prior.Uniform(minimum=0, maximum=0.8, name=None, 
        #                                                                     latex_label=None, unit=None, boundary=None), 
        #                                     z_prior=bilby.core.prior.Uniform(minimum=-1, maximum=1, name=None, 
        #                                                                     latex_label=None, unit=None, boundary=None), 
        #                                     name='chi_2', latex_label='$\\chi_2$', unit=None, boundary=None, 
        #                                     minimum=-0.8, maximum=0.8)

        priors["chi_1"] = bilby.core.prior.Uniform(
    minimum=-0.8, maximum=0.8, name="chi_1", latex_label=r"$\chi_1$"
    )
        priors["chi_2"] = bilby.core.prior.Uniform(
            minimum=-0.8, maximum=0.8, name="chi_2", latex_label=r"$\chi_2$"
        )
        inj_params = dicc.get("inyection", {})
        priors['theta_jn'] = float(inj_params.get('theta_jn', 0.9))
        priors['ra'] = float(inj_params.get('ra', 1.375))
        priors['dec'] = float(inj_params.get('dec', -1.21))
        priors['psi'] = float(inj_params.get('psi', 0.5))
        priors['phase'] = bilby.core.prior.Uniform(
            minimum=0.0,
            maximum=2 * np.pi,
            name='phase',
            latex_label=r'$\phi$',
            boundary='periodic'
        )
        priors['luminosity_distance'] = float(inj_params.get('luminosity_distance', dist_inj))
        priors['geocent_time'] = bilby.core.prior.Uniform(
            minimum=tc_inj - 0.05,
            maximum=tc_inj + 0.05,
            name='geocent_time',
            latex_label=r'$t_c$'
        )
        targ_keys = [
            "chirp_mass",
            "mass_ratio",
            "chi_1",
            "chi_2",
            "luminosity_distance",
            "geocent_time",
            "phase",
            "theta_jn",
            "ra",
            "dec",
            "psi"
        ]
    else:
        # Precessing model
        priors = bilby.gw.prior.BBHPriorDict(aligned_spin=False)

        priors['chirp_mass'] = bilby.gw.prior.UniformInComponentsChirpMass(
            minimum=10,
            maximum=100,
            name='chirp_mass',
            latex_label=r'$\mathcal{M}$'
        )
        #         priors['chirp_mass'] = bilby.gw.prior.UniformInComponentsChirpMass(
        #     minimum=max(5.0, mc_inj * 0.4),
        #     maximum=mc_inj * 1.8,
        #     name='chirp_mass',
        #     latex_label=r'$\mathcal{M}$'
        # )
        priors['mass_ratio'] = bilby.gw.prior.UniformInComponentsMassRatio(
            minimum=1.0 / 8.0,
            maximum=1.0,
            name='mass_ratio',
            latex_label=r'$q$'
        )
        priors['a_1'] = bilby.core.prior.Uniform(0.0, 0.8, name='a_1', latex_label=r'$a_1$')
        priors['a_2'] = bilby.core.prior.Uniform(0.0, 0.8, name='a_2', latex_label=r'$a_2$')
        priors['tilt_1'] = bilby.core.prior.Sine(name='tilt_1', latex_label=r'$\theta_1$')
        priors['tilt_2'] = bilby.core.prior.Sine(name='tilt_2', latex_label=r'$\theta_2$')
        priors['phi_12'] = bilby.core.prior.Uniform(0.0, 2 * np.pi, name='phi_12', latex_label=r'$\Delta\phi$', boundary='periodic')
        priors['phi_jl'] = bilby.core.prior.Uniform(0.0, 2 * np.pi, name='phi_jl', latex_label=r'$\phi_{JL}$', boundary='periodic')
        priors['theta_jn'] = 0.0
        priors['ra'] = 0.0
        priors['dec'] = 0.0
        priors['psi'] = 0.0
        priors['phase'] = bilby.core.prior.Uniform(
            minimum=0.0,
            maximum=2 * np.pi,
            name='phase',
            latex_label=r'$\phi$',
            boundary='periodic'
        )
        priors['luminosity_distance'] = dist_inj
        priors['geocent_time'] = bilby.core.prior.Uniform(
            minimum=tc_inj - 0.05,
            maximum=tc_inj + 0.05,
            name='geocent_time',
            latex_label=r'$t_c$'
        )
        targ_keys = [
            "chirp_mass",
            "mass_ratio",
            "a_1",
            "a_2",
            "tilt_1",
            "tilt_2",
            "phi_12",
            "phi_jl",
            "luminosity_distance",
            "geocent_time",
            "phase"
        ]

    print(' ')
    print('*' * 50)
    print('Priors loaded successfully')
    print('*' * 50)
    print(' ')
    return targ_keys, priors

