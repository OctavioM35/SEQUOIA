import h5py
import numpy as np

# ---------------------------------------------------------
# 4. Check chirp mass
# ---------------------------------------------------------

def check_chirp_mass(zenodo_file):
    """
    Check whether the median chirp mass of an event is high enough for the
    surrogate models.

    The function opens the LVK Zenodo file, takes the first analysis in the
    file and computes the median of its chirp mass posterior samples. The
    event is accepted only if this median is above 25 solar masses, since
    lighter systems have long signals that lead to excessive truncation
    with the surrogates (DANSur and NRSur7dq4).

    INPUTS:
    -------

    -zenodo_file: Path to the LVK Zenodo .h5 file of the event.

    OUTPUTS:
    --------

    -True if the median chirp mass is greater than 25 solar masses.

    -False otherwise.
    """
    with h5py.File(zenodo_file, "r") as f:
        ap_key = list(f.keys())[0]
        chirp = np.median(f[ap_key]["posterior_samples"]["chirp_mass"])
        if chirp <= 25 :
            return False
        else:
            return True