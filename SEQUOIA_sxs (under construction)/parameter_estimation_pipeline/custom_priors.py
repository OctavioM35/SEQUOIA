import h5py
import bilby
from bilby.gw.prior import Uniform, Constraint, UniformInComponentsChirpMass, UniformInComponentsMassRatio
import numpy as np


def load_custom_priors(dicc):
    
            """
            Establish custom priors and targ keys for the PE

            INPUTS:
            ------

            -dicc: The created dicctionay in config.py 

            Returns
            ------

            -targ_keys: Names of the source parameters used by the waveform generator.

            -priors: Prior distributions for the parameters sampled during the Bayesian
            inference.
            """
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
                            priors['theta_jn'] = 1.
                            priors['ra'] = 0.
                            priors['dec'] =0.
                            priors['psi'] =0.
                
                            # priors['mass_ratio'] = bilby.gw.prior.UniformInComponentsMassRatio(minimum=1/8, maximum=1.0, name='mass_ratio', latex_label='$q$', unit=None, boundary=None, equal_mass=False)

                            targ_keys = [
                                "chirp_mass",
                                "mass_ratio",
                                "chi_1",
                                "chi_2",
                                "luminosity_distance",
                                "geocent_time"
                            ]                
            print(' ')
            print('*' *50)
            print('Priors loaded')
            print('*' *50)
            print(' ')
            return targ_keys, priors

