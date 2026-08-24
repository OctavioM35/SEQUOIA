import h5py
import bilby
from bilby.gw.prior import Uniform, Constraint, UniformInComponentsChirpMass, UniformInComponentsMassRatio
import numpy as np


def load_custom_priors():

    priors = bilby.gw.prior.BBHPriorDict(aligned_spin=True)

    
    f'''
    Write custom priors


    priors['luminosity_distance']=bilby.gw.prior.UniformSourceFrame(minimum=100.0, maximum=13000, 
        cosmology='Planck15', name='luminosity_distance', latex_label='$d_L$', unit='Mpc', 
        boundary=None)
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
    '''
    targ_keys = bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters(priors.sample())[0]
    targ_keys = set(targ_keys.keys())

    return targ_keys,priors