import os
from config import *
from parameter_estimation_pipeline.custom_priors import load_custom_priors
from parameter_estimation_pipeline.run_inference_bilby import run_inference
# from parameter_estimation_pipeline.plot_manual_corner import  plot_manual_corner
from utils.evaluate_selected_surrogate import evaluate_selected_surrogate
from parameter_estimation_pipeline.mismatch import *
from utils.synthetic_waveform import synthetic_waveform
import shutil

def process_event(sim,surrogate_model, results):

    try:    
        event_dir = os.path.join(sim)
        # ---------------------------------------------------------
        # 1. Create output directory
        # ---------------------------------------------------------

        outdir = os.path.join(results,f"outdir_NN_{sim}")
        os.makedirs(outdir, exist_ok=True)
        precession_or_not = evaluate_selected_surrogate(surrogate_model)

        # ---------------------------------------------------------
        # 2. Skip event if corner comparison already exists
        # ---------------------------------------------------------
        if skip:

            existing_corner = any(
                "corner_comparison" in filename
                for filename in os.listdir(outdir)
            )


            if existing_corner:

                print(
                    f"Skipping {sim}: "
                    "corner_comparison already exists."
                )

                return None

      

        # ---------------------------------------------------------
        # 5. Load LVK configuration
        # ---------------------------------------------------------

        dicc["precession_model"] = precession_or_not

        targ_keys, priors = load_custom_priors(dicc)
  
        # ---------------------------------------------------------
        # 6. Check priors
        # ---------------------------------------------------------


        if targ_keys is None or priors is None:

            print(
                f"Priors are not valid for {sim}. Skipping event..."
            )
            if os.path.exists(outdir):
                    shutil.rmtree(outdir)
            return "Priors not valid."


        # ---------------------------------------------------------
        # 7. Load interferometers
        # ---------------------------------------------------------

        dicc["sxs"] = sim
        ifos = synthetic_waveform(dicc)


        if ifos is None:

            print(
                f"This signal does not exist "
                f"for event {sim}."
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
            sim
        )

        # ---------------------------------------------------------
        # 9. Plot LVK and obtained posterior distributions
        # ---------------------------------------------------------

        # print("=" * 60)
        # print("Plotting results...")
        # print("=" * 60)
        # print('Calculating mismatch...')
        # O = overlap(dicc["signal"][1], dicc["generated_waveform"],  dicc["sampling-frequency"])
        # loss = L2(dicc["signal"][1], dicc["generated_waveform"],  dicc["sampling-frequency"])

        # print(f"Overlap = {O:.8f}")
        # print(f"L2      = {loss:.8f}")

        print(f"Event {sim} processed successfully.")

        
    
        return None
    
    except ZeroDivisionError as e:


        if (
            isinstance(e, ValueError)
            and str(e) == 'f_ref cannot be lower than f_low.'
        ):
            problem = 'This surrogate does not allow f_ref < f_low.'
            print(f"Exception while processing {sim}: DANSur does not allow f_ref < f_low.")

        else:
            print(f"Exception while processing {sim}: {e}")
            problem = str(e)

        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        return problem