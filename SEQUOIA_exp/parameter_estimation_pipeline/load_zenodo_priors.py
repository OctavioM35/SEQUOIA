import bilby
import h5py
from bilby.gw.prior import *
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine, Cosine
import ast



def load_zenodo_priors(zenodo_file):
    
    with h5py.File(zenodo_file, "r") as f:
            key0 = list(f.keys())[0]
            
            print('Loading LVK  priors')
            available_priors = f[key0]["priors"]["analytic"]


            priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)
            # priors.clear()
            # priors.conversion_function = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters
            # priors['mass_ratio'] = eval(available_priors['mass_ratio'][:][0])
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
            # PROBAR CON LOS 2 DE ABAJO PUESTOS (O AL MENOS SIN LA MASA CHIRP)
            priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
            priors['chirp_mass'] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )
            priors['mass_ratio'] = eval(f[key0]["priors"]["analytic"]["mass_ratio"][:][0] ) 
            # priors['mass_1'] = eval(f[key0]["priors"]["analytic"]["mass_1"][:][0] )
            # priors['mass_2'] = eval(f[key0]["priors"]["analytic"]["mass_2"][:][0] ) 
            # priors.conversion_function = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters

            # if approximant == True:  
            targ_keys = [
                    "chi_1",
                    "chi_2",
                    "luminosity_distance",
                    "geocent_time",
                    "chirp_mass",
                    "mass_ratio"
                ]
            # try:
            #     priors["dec"] = eval(f[key0]["priors"]["analytic"]["dec"][:][0])
            #     priors["ra"] = eval(f[key0]["priors"]["analytic"]["ra"][:][0])
            # except KeyError:
            #     pass
            #     targ_keys = [
            #         "mass_1",
            #         "mass_2",
            #         "a_1",
            #         "a_2",
            #         "tilt_1",
            #         "tilt_2",
            #         "phi_12",
            #         "phi_jl",
            #         "luminosity_distance",
            #         "theta_jn",
            #         "psi",
            #         "phase",
            #         "geocent_time",
            #         "ra",
            #         "dec"
            #     ]

            # ==========================================================
            # TARG_KEYS
            # ==========================================================



            # print(targ_keys)
            # raise Exception
            # print('priors: ')
            # prior = priors["chirp_mass"]

            # array = prior.sample(1000)

            # print(array)
            # raise Exception
            targ_keys = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(priors.sample())[0]
            targ_keys = set(targ_keys.keys())
    print(' ')
    print('*' *50)
    print('Priors loaded')
    print('*' *50)
    print(' ')
    return priors, targ_keys

