import os
import shutil

from config import *
from sequoiasxs.priors.custom_priors import load_custom_priors
from sequoiasxs.pipeline.run_inference_bilby  import run_inference
from utils.evaluate_selected_surrogate import evaluate_selected_surrogate
from utils.synthetic_waveform import synthetic_waveform


def process_event(sim,surrogate_model, results):
    """
    Run the full parameter-estimation workflow for a single SXS simulation.

    The function creates the output directory of the simulation, loads the
    custom priors, builds the synthetic H1 and L1 data from the SXS waveform
    and runs the Bayesian inference with the selected surrogate model. If the
    simulation was already processed (skip option enabled and a results file
    already present in the output directory), it is skipped. If anything fails
    (invalid priors, missing data, or an exception during the inference), the
    output directory is removed and a string describing the problem is returned.

    INPUTS:
    -------

    -sim: Identifier of the SXS simulation (e.g. "SXS:BBH:1156").

    -surrogate_model: Name of the surrogate model used for the inference.
    It is also used to determine whether the model includes precession.

    -results: Path to the directory where the output folder of the
    simulation (outdir_NN_<sim>) is created.

    OUTPUTS:
    --------

    -None if the simulation was processed successfully or skipped.

    -A string describing the problem if it could not be processed
    (e.g. "Priors not valid.", "Interferometer data missing.", or the
    exception message).
    """

    try:    
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
                "Obtained vs inyected parameters for surrogate" in filename
                for filename in os.listdir(outdir)
            )


            if existing_corner:

                print(
                    f"Skipping {sim}: "
                    "Event already exists."
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
    
    except Exception as e:


        if (
            isinstance(e, ValueError)
            and str(e) == 'f_ref cannot be lower than f_low.'
        ):
            problem = 'This surrogate does not allow f_ref < f_low.'
            print(f"Exception while processing {sim}: This surrogate does not allow f_ref < f_low.")

        else:
            print(f"Exception while processing {sim}: {e}")
            problem = str(e)

        if os.path.exists(outdir):
            shutil.rmtree(outdir)

        return problem