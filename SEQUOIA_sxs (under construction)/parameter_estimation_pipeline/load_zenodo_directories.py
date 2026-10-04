import os

# ---------------------------------------------------------
# 3. Load Zenodo data
# ---------------------------------------------------------

def load_zenodo(data_folder,folder,event_dir):
        zenodo_file = None
        for f in os.listdir(event_dir):
            if f.endswith(folder + ".h5") or f.endswith(folder +".hdf5") or f.endswith("cosmo.h5") or f.endswith("cosmo.hdf5"):
                zenodo_file = os.path.join(event_dir, f)

                break

        if zenodo_file is None:
            return None
        zenodo_file = os.path.join(data_folder,folder, zenodo_file) 
        return zenodo_file 