#main.py

import os
import sys
import sxs

from config import *
from pathlib import Path
from sequoiasxs.pipeline.process_event import process_event
from sequoiasxs.supported_surrogates.install_surrogates import install_surrogates


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def main():
    """
    Main program for SEQUOIASXS
    """
    events_with_issues = []

    if install_surrogates == True:
           install_surrogates()
  

    if run_particular_event:

        sims = ['SXS:' + sxs_simulation]
    else:

        sims = [ ]
        for sxs_id in [
    'SXS:BBH:0110v2.0' , 'SXS:BBH:0237v2.0' , 'SXS:BBH:0926','SXS:BBH:1156','SXS:BBH:1355','SXS:BBH:1419','SXS:BBH:1509']: 
                sims.append(sxs_id)
    total_events = len(sims)
    print(f"Found {total_events} event(s)")
    for i, sim in enumerate(sims, start=1):

            print("\n" + "=" * 60)
            print(f"Processing event {sim}")
            print(f"Event {i}/{total_events}")
            print("=" * 60)

                    
            problem = process_event(sim,surrogate_model,results)

            if problem is not None:
                    events_with_issues.append(
                        (sim, problem)
                    )


    print("\n" + "=" * 60)
    print("Pipeline execution completed")
    print("=" * 60)

    print("Events with issues: ")

    for folder, reason in events_with_issues:
                print(f"  - {folder}: {reason}")
    print('Number of events with issues: ' , len(events_with_issues))


if __name__ == "__main__":
    main()