import bilby
import h5py
from bilby.gw.prior import *
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine, Cosine
import ast


# ---------------------------------------------------------
# 5.1 Load priors
# ---------------------------------------------------------


def load_zenodo_priors(zenodo_file,dicc):
    
    with h5py.File(zenodo_file, "r") as f:
            
            key0 = list(f.keys())[0]
            print('Loading LVK  priors')
            
            if dicc["precession_model"] == False:
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
                            priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
                            priors['chirp_mass'] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )
                            priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=0.125, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)
                            priors['mass_1'] = eval(f[key0]["priors"]["analytic"]["mass_1"][:][0] )
                            priors['mass_2'] = eval(f[key0]["priors"]["analytic"]["mass_2"][:][0] ) 
                    

                            try:
                                priors['ra'] = eval(f[key0]["priors"]["analytic"]["ra"][:][0] ) 
                                priors['dec'] = eval(f[key0]["priors"]["analytic"]["dec"][:][0] ) 
                            except Exception:
                                   pass



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

            else:
                            priors = bilby.gw.prior.BBHPriorDict(aligned_spin=False)
                            priors['a_1'] = bilby.core.prior.Uniform(0, 0.8, name='a_1')
                            priors['a_2'] = bilby.core.prior.Uniform(0, 0.8, name='a_2')
                
                
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
                            priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
                            priors['chirp_mass'] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )
                            priors["mass_ratio"] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=1/4, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)
                            priors['mass_1'] = eval(f[key0]["priors"]["analytic"]["mass_1"][:][0] )
                            priors['mass_2'] = eval(f[key0]["priors"]["analytic"]["mass_2"][:][0] ) 
                            
                            #Si se incluye precesión
                            priors['tilt_1'] = eval(f[key0]["priors"]["analytic"]["tilt_1"][:][0] )
                            priors['tilt_2'] = eval(f[key0]["priors"]["analytic"]["tilt_2"][:][0] ) 
                            priors['phi_12'] = eval(f[key0]["priors"]["analytic"]["phi_12"][:][0] )
                            priors['phi_jl'] = eval(f[key0]["priors"]["analytic"]["phi_jl"][:][0] ) 



                            targ_keys = [
                                "chirp_mass",
                                "mass_ratio",
                                "a_1",
                                "a_2",
                                "luminosity_distance",
                                "theta_jn",
                                "phase",
                                "ra",
                                "dec",
                                "geocent_time",
                                "psi",
                                #Si se incluye precesión:
                                "tilt_1",
                                "tilt_2",
                                "phi_12",
                                "phi_jl"
                            ]


                          
                          
    print(' ')
    print('*' *50)
    print('Priors loaded')
    print('*' *50)
    print(' ')
    return priors, targ_keys

