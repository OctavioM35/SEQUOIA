import os

# ---------------------------------------------------------
# 3. Load Zenodo data
# ---------------------------------------------------------

def load_zenodo(data_folder,folder,event_dir):
        """
        Locate the LVK Zenodo file of an event inside its directory.

        The function scans the event directory and selects the first file whose
        name ends with "<folder>.h5" or "<folder>.hdf5" (the event file) or with
        "cosmo.h5" or "cosmo.hdf5" (the version of the file that includes the
        cosmology). The scan stops at the first match. If no file matches, None
        is returned.

        INPUTS:
        -------

        -data_folder: Path to the root directory containing the event
        directories.

        -folder: Name of the event directory (e.g. "GW230825_041334"). It is
        also used as the suffix of the expected file name.

        -event_dir: Path to the event directory that is scanned.

        OUTPUTS:
        --------

        -Path to the Zenodo file, built as data_folder/folder/<path>, if it is
        found.

        -None if no matching file exists in the directory.
        """        
        zenodo_file = None
        for f in os.listdir(event_dir):
            if f.endswith(folder + ".h5") or f.endswith(folder +".hdf5") or f.endswith("cosmo.h5") or f.endswith("cosmo.hdf5"):
                zenodo_file = os.path.join(event_dir, f)

                break

        if zenodo_file is None:
            return None
        zenodo_file = os.path.join(data_folder,folder, zenodo_file) 
        return zenodo_file 