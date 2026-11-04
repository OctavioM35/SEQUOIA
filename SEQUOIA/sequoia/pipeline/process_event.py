import os
import shutil

#Load utils
from config import *
from sequoia.utils.chirp_mass import check_chirp_mass 
from sequoia.utils.generate_dictionary import load_zenodo_configuration


#Load data
from sequoia.data.load_zenodo_directories import load_zenodo
from sequoia.data.load_interferometer import load_ifos

#Load priors
from sequoia.priors.custom_priors import load_custom_priors
from sequoia.priors.load_zenodo_priors import load_zenodo_priors

#Load surrogate
from sequoia.supported_surrogates.evaluate_selected_surrogate import evaluate_selected_surrogate

#Load inference and plots
from sequoia.pipeline.run_inference_bilby import run_inference
from sequoia.plotting.plot_manual_corner import  plot_manual_corner



def process_event(folder,surrogate_model, results):
    """
    Run the full parameter-estimation workflow for a single LVK event.

    The function creates the output directory of the event, loads its LVK
    Zenodo file, checks that the chirp mass is high enough for the surrogate,
    and loads the configuration and priors (from the Zenodo file if the
    automatic option is enabled, or from the custom priors otherwise). It
    then loads the H1, L1 and V1 data, runs the Bayesian inference with the
    selected surrogate model, and produces the corner plots, including the
    comparison with the LVK posteriors. If the event was already processed
    (skip option enabled and a corner comparison file present in the output
    directory), it is skipped. If anything fails (missing Zenodo file, chirp
    mass too low, invalid priors, missing interferometer data, or an
    exception during the inference), the output directory is removed and a
    string describing the problem is returned.

    INPUTS:
    -------

    -folder: Name of the event directory (e.g. "GW230825_041334"), located
    inside data_folder.

    -surrogate_model: Name of the surrogate model used for the inference.
    It is also used to determine whether the model includes precession and
    to set the label of the outputs.

    -results: Path to the directory where the output folder of the event
    (outdir_NN_<folder>) is created.

    OUTPUTS:
    --------

    -None if the event was processed successfully or skipped.

    -A string describing the problem if it could not be processed
    (e.g. "Zenodo file missing.", "Chirp mass too low: will lead to
    excessive truncation.", "Priors not valid.", "Interferometer data
    missing.", or the exception message).
    """
    try:    
        event_dir = os.path.join(data_folder, folder)
        # ---------------------------------------------------------
        # 1. Create output directory
        # ---------------------------------------------------------
        outdir = os.path.join(
                results,
                f"outdir_NN_{folder}"
            )


        outdir = os.path.join(results,f"outdir_NN_{folder}")
        os.makedirs(outdir, exist_ok=True)
        precession_or_not, label = evaluate_selected_surrogate(surrogate_model)

        # ---------------------------------------------------------
        # 2. Skip event if corner comparison already exists
        # ---------------------------------------------------------
        print(outdir)
        if skip:

            existing_corner = any(
                "corner_comparison" in filename
                for filename in os.listdir(outdir)
            )


            if existing_corner:

                print(
                    f"Skipping {folder}: "
                    "corner_comparison already exists."
                )

                return None

        # ---------------------------------------------------------
        # 3. Load Zenodo data
        # ---------------------------------------------------------
        zenodo_file = load_zenodo(
            data_folder,
            folder,
            event_dir
        )

        if zenodo_file is None:

            print(
                f"There is no Zenodo data "
                f"for event {folder}. Skipping event..."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "Zenodo file missing."

        # ---------------------------------------------------------
        # 4. Check chirp mass
        # ---------------------------------------------------------

        if not check_chirp_mass(zenodo_file):

            print(
                "Chirp mass too low -will lead to "
                "excesive truncation in Dansur and NRSur-. Skipping event..."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "Chirp mass too low: will lead to excesive truncation."

        # ---------------------------------------------------------
        # 5. Load LVK configuration
        # ---------------------------------------------------------

        if automatic:

            dicc = load_zenodo_configuration(
                zenodo_file
            )
            dicc["precession_model"] = precession_or_not
            dicc["surrogate_model"] = surrogate_model
            dicc["label"] = label
        # ---------------------------------------------------------
        # 5.1 Load priors
        # ---------------------------------------------------------


            priors, targ_keys = load_zenodo_priors(zenodo_file,dicc)
            
        else:
            dicc["precession_model"] = precession_or_not

            targ_keys, priors = load_custom_priors()

   
        # ---------------------------------------------------------
        # 6. Check priors
        # ---------------------------------------------------------


        if targ_keys is None or priors is None:

            print(
                f"Priors are not valid for {folder}. Skipping event..."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "Priors not valid."


        # ---------------------------------------------------------
        # 7. Load interferometers
        # ---------------------------------------------------------

        
        ifos = load_ifos(
                event_dir,
                folder,
                dicc,
                zenodo_file,
            )




        if ifos is None:

            print(
                f"There is no interferometer data "
                f"for event {folder}."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "Interferometer data missing."

        # ---------------------------------------------------------
        # 8. Run inference
        # ---------------------------------------------------------

        run_inference(
            ifos,
            dicc,
            targ_keys,
            priors,
            outdir,
            folder
        )

        # ---------------------------------------------------------
        # 9. Plot LVK and obtained posterior distributions
        # ---------------------------------------------------------

        print("=" * 60)
        print("Plotting results...")
        print("=" * 60)

        plot_manual_corner(
            folder,
            zenodo_file,
            event_dir,
            outdir,
            dicc
        )

        print(f"Event {folder} processed successfully.")
        return None
    
    except Exception as e:


        if (
            isinstance(e, ValueError)
            and str(e) == 'f_ref cannot be lower than f_low.'
        ):
            problem = 'This surrogate does not allow f_ref < f_low.'
            print(f"Exception while processing {folder}: DANSur does not allow f_ref < f_low.")

        else:
            print(f"Exception while processing {folder}: {e}")
            problem = str(e)

        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        return problem