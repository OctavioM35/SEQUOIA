import bilby
import h5py
import ast

from bilby.gw.prior import *
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine, Cosine

# ---------------------------------------------------------
# 5.1 Load priors
# ---------------------------------------------------------


def load_zenodo_priors(zenodo_file,dicc):
    """
    Load the prior distributions used by the LVK for an event from its
    Zenodo file.

    The function reads from the "priors/analytic" group of the first
    analysis in the Zenodo file the priors of the geocentric time (or, if
    it does not exist, of the L1 or H1 time), luminosity distance,
    theta_jn, polarization angle, phase, chirp mass, component masses and,
    when available, right ascension and declination. These priors are
    stored in the file as strings and are evaluated into Bilby prior
    objects, so no distribution is built from scratch for them.

    Two priors are not taken from the file. The mass-ratio prior is
    uniform in the components, with a lower limit that depends on the
    validity domain of the selected surrogate (1/8 for DANSurHM and
    NRHybSur3dq8, 1/4 for NRSur7dq4). The aligned-spin priors (chi_1,
    chi_2) are set to the Bilby aligned-spin prior with a spin magnitude
    uniform in [-0.8, 0.8] and a uniform z-component in [-1, 1], starting
    from the default BBHPriorDict with aligned spins.

    INPUTS:
    -------

    -zenodo_file: Path to the LVK Zenodo .h5 file of the event.

    -dicc: Configuration dictionary. The key "label" (name of the
    surrogate model) is used to select the mass-ratio prior.

    OUTPUTS:
    --------

    -priors: bilby.gw.prior.BBHPriorDict with the priors of the inference.

    -targ_keys: List with the names of the sampled parameters: chirp_mass,
    mass_ratio, chi_1, chi_2, luminosity_distance, theta_jn, phase, psi,
    geocent_time, ra and dec.

    -(None, None) if no time prior (geocent_time, L1_time or H1_time) is
    found in the file.
    """    
    with h5py.File(zenodo_file, "r") as f:
            
                            key0 = list(f.keys())[0]
                            print('Loading LVK  priors')
            
                            #Sin precesión                        
                            priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)     
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
                            # # # # PROBAR CON LOS 2 DE ABAJO PUESTOS (O AL MENOS SIN LA MASA CHIRP)
                            priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
                            priors['chirp_mass'] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )
                            if dicc["label"] == 'DANSurHM':
                                priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=1/8, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)
                            elif dicc["label"] =='NRHybSur3dq8':
                                priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=1/8, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)
                            elif dicc["label"] =='NRSur7dq4':
                                priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=1/4, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)
                            priors['mass_1'] = eval(f[key0]["priors"]["analytic"]["mass_1"][:][0] )
                            priors['mass_2'] = eval(f[key0]["priors"]["analytic"]["mass_2"][:][0] ) 
                            # priors['tilt_1'] = eval(f[key0]["priors"]["analytic"]["tilt_1"][:][0] )
                            # priors['tilt_2'] = eval(f[key0]["priors"]["analytic"]["tilt_2"][:][0] ) 
                            # priors['phi_12'] = eval(f[key0]["priors"]["analytic"]["phi_12"][:][0] )
                            # priors['phi_jl'] = eval(f[key0]["priors"]["analytic"]["phi_jl"][:][0] ) 
                            # priors['azimuth'] = eval(f[key0]["priors"]["analytic"]["azimuth"][:][0] )
                            # priors['zenith'] = eval(f[key0]["priors"]["analytic"]["zenith"][:][0] ) 
                            # priors['time_jitter'] = eval(f[key0]["priors"]["analytic"]["time_jitter"][:][0] ) 

                            try:
                                priors['ra'] = eval(f[key0]["priors"]["analytic"]["ra"][:][0] ) 
                                priors['dec'] = eval(f[key0]["priors"]["analytic"]["dec"][:][0] ) 
                            except Exception:
                                   pass

                            # analytic_priors = f[key0]["priors"]["analytic"]

                            # for key in analytic_priors.keys():

                            #     if "recalib" in key.lower():

                            #         priors[key] = eval(
                            #             analytic_priors[key][:][0]
                            #         )

                            #         print(
                            #             f"Loaded calibration prior {key}: "
                            #             f"{priors[key]}"
                            #         )


                            targ_keys = [
                                "chirp_mass",
                                "mass_ratio",
                                "chi_1",
                                "chi_2",
                                "luminosity_distance",
                                "theta_jn",
                                "phase",
                                "psi",
                                "geocent_time",
                                "ra",
                                "dec",
                            ]
                            # targ_keys = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(priors.sample())[0]

                            # targ_keys = set(targ_keys.keys())
                          
                          
    print(' ')
    print('*' *50)
    print('Priors loaded')
    print('*' *50)
    print(' ')
    return priors, targ_keys

