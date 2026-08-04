import bilby
import h5py
from bilby.gw.prior import *
from astropy.cosmology import LambdaCDM
from bilby.core.prior import Sine
import ast

def rewrite(s,parameter,envelope):
    #  { H1:20,L1:20, waveform: 11.0 }
    d={}
    config = s.strip('{}  ,').split(',')

    # if parameter == 'spline-calibration-envelope-dict':
    #     for i in config:
    #         param=i.split(':')
    #         d[param[0]] = envelope[param[0]]
        
    #     return d

    # else:
    for i in config:
            param=i.split(':')
            try:
                d[param[0]] = ast.literal_eval(param[1])
            except (ValueError, SyntaxError):
                d[param[0]] = param[1]
            except IndexError:
                return s    
    return d

def load_zenodo_configuration(zenodo_file):
    
    with h5py.File(zenodo_file, "r") as f:
            key0 = list(f.keys())[0]
            print('=' *50)
            print('Loading Zenodo configuration...')
            print(' ')
            config = f[key0]["config_file"]["config"]
            dicc = {}
            envelope = f[key0]["calibration_envelope"]
            for parameter, ds in config.items():
                s = ds[0].decode("utf-8").strip()

                try:
                    dicc[parameter] = ast.literal_eval(s)
                except (ValueError, SyntaxError):
                    if s.startswith('{') and s.endswith('}'):
                        dicc[parameter] = rewrite(s,parameter,envelope)
                    else:
                        dicc[parameter] = s
            print('Loading completed.')

    
            print('Loading priors')
            available_priors = f[key0]["priors"]["analytic"]


            priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)
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
            # print(type(priors))
            # raise Exception
            # for prior in list(available_priors):
            #         if prior =='H1_time':
            #             priors["geocent_time"] = eval(f[key0]["priors"]["analytic"][prior][:])
            #         elif prior =='mass_1' or prior=='mass_2':
            #             continue
            #         else:
            #             priors[prior] = eval(f[key0]["priors"]["analytic"][prior][:])


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
            # priors["a_1"] = eval(f[key0]["priors"]["analytic"]["a_1"][:][0] )
            # priors["a_2"] = eval(f[key0]["priors"]["analytic"]["a_2"][:][0] )
            # priors["chirp_mass"] = eval(f[key0]["priors"]["analytic"]["chirp_mass"][:][0] )
            priors["luminosity_distance"] = eval(f[key0]["priors"]["analytic"]["luminosity_distance"][:][0] )
            priors["mass_1"] = eval(f[key0]["priors"]["analytic"]["mass_1"][:][0] )
            priors["mass_2"] = eval(f[key0]["priors"]["analytic"]["mass_2"][:][0] )
            # priors["mass_ratio"] = eval(f[key0]["priors"]["analytic"]["mass_ratio"][:][0] )
            # priors["phase"] = eval(f[key0]["priors"]["analytic"]["phase"][:][0] )
            # priors["phi_12"] = eval(f[key0]["priors"]["analytic"]["phi_12"][:][0] )
            # priors["phi_jl"] = eval(f[key0]["priors"]["analytic"]["phi_jl"][:][0] )
            priors["psi"] = eval(f[key0]["priors"]["analytic"]["psi"][:][0] )
            priors["theta_jn"] = eval(f[key0]["priors"]["analytic"]["theta_jn"][:][0] )
            # priors["tilt_1"] = eval(f[key0]["priors"]["analytic"]["tilt_1"][:][0] )
            # priors["tilt_2"] = eval(f[key0]["priors"]["analytic"]["tilt_2"][:][0] )
            # # priors["recalib_H1_amplitude_0"] = eval(f[key0]["priors"]["analytic"]["recalib_H1_amplitude_0"][:][0] )
            
            try:
                   
                   priors["dec"] = eval(f[key0]["priors"]["analytic"]["dec"][:][0])
                   priors["ra"] = eval(f[key0]["priors"]["analytic"]["ra"][:][0])
            except Exception:
                    pass

                  
            targ_keys = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(priors.sample())[0]
            targ_keys = set(targ_keys.keys())
            # print(targ_keys)
            # raise Exception
            # print('priors: ')
            # prior = priors["chirp_mass"]

            # array = prior.sample(1000)

            # print(array)
            # raise Exception
    print(' ')
    print('*' *50)
    print('Priors loaded')
    print('*' *50)
    print(' ')
    return dicc, priors, targ_keys

