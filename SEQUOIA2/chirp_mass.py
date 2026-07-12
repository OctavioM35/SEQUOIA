import h5py
import numpy as np

def chirp_mass(zenodo_file):
    with h5py.File(zenodo_file, "r") as f:
        ap_key = list(f.keys())[0]
        chirp = np.median(f[ap_key]["posterior_samples"]["chirp_mass"])
        print(chirp)
        if chirp < 22 or chirp > 100 :
            return False
        else:
            return True