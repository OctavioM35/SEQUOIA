import os
from config import *
from parameter_estimation_pipeline.load_zenodo_directories import load_zenodo
from parameter_estimation_pipeline.custom_priors import load_custom_priors
from parameter_estimation_pipeline.load_interferometer import load_ifos
from parameter_estimation_pipeline.run_inference_bilby import run_inference
from parameter_estimation_pipeline.plot_manual_corner import  plot_manual_corner
from utils.chirp_mass import chirp_mass 
from utils.generate_dictionary import load_zenodo_configuration
from utils.surrogates import surrogates

from parameter_estimation_pipeline.load_zenodo_priors import load_zenodo_priors
from utils.synthetic_waveform import synthetic_waveform
import shutil

def process_event(folder):
    # if surrogate_model not in surrogates:
    #      raise KeyError('Surrogate not supported by')
    # try:    
        event_dir = os.path.join(data_folder, folder)

        # ---------------------------------------------------------
        # 1. Create output directory
        # ---------------------------------------------------------
        if generate_synthetic_signal == False:
            outdir = os.path.join(
                results,
                f"outdir_NN_{folder}"
            )
        elif generate_synthetic_signal == True:
            outdir = os.path.join(
                results,
                f"outdir_synthetic_NN_{folder}"
            )
        else:
            raise KeyError('generate_synthetic_signal must either be True or false.')
        os.makedirs(outdir, exist_ok=True)

        # ---------------------------------------------------------
        # 2. Skip event if corner comparison already exists
        # ---------------------------------------------------------

        if skip:

            existing_corner = any(
                filename.startswith("corner_comparison")
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

        if not chirp_mass(zenodo_file) and approximant == False:

            print(
                "DANSur has not been trained for handling "
                "this particular event "
                "(low or high mass event). Skipping event..."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "DANSur does not support this event (low or high mass event)."

        # ---------------------------------------------------------
        # 5. Load LVK configuration
        # ---------------------------------------------------------


        if automatic:

            dicc = load_zenodo_configuration(
                zenodo_file
            )


        # ---------------------------------------------------------
        # 5.1 Load priors
        # ---------------------------------------------------------


            priors, targ_keys = load_zenodo_priors(zenodo_file,dicc)
            
        else:

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


        if generate_synthetic_signal == False:
            ifos = load_ifos(
                event_dir,
                folder,
                dicc,
                zenodo_file,
            )
        else:
            ifos = synthetic_waveform(dicc,zenodo_file)




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
    
    # except Exception as e:


    #     if (
    #         isinstance(e, ValueError)
    #         and str(e) == 'f_ref cannot be lower than f_low.'
    #     ):
    #         problem = 'DANSur does not allow f_ref < f_low.'
    #         print(f"Exception while processing {folder}: DANSur does not allow f_ref < f_low.")

    #     else:
    #         print(f"Exception while processing {folder}: {e}")
    #         problem = str(e)

    #     if os.path.exists(outdir):
    #         shutil.rmtree(outdir)

        return problem