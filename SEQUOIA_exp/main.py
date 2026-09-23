#main.py



import os
from config import *
from process_event import process_event
import sys
from pathlib import Path
from utils.install_surrogates import install_surrogates

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def main():

    if install_surrogates is True:
          install_surrogates()

    events_with_issues = []
    
    folders = [
        folder
        for folder in os.listdir(data_folder)
        if os.path.isdir(os.path.join(data_folder, folder))
    ]

    if run_particular_event:

        folders = [particular_event]

    else:

        folders = [
            folder
            for folder in os.listdir(data_folder)
            if os.path.isdir(os.path.join(data_folder, folder))
        ]

    total_events = len(folders)

    print(f"Found {total_events} event(s)")
   
    for i, folder in enumerate(folders, start=1):
            # if folder in ruido:
            #     continue
            print("\n" + "=" * 60)
            print(f"Processing event {folder}")
            print(f"Event {i}/{total_events}")
            print("=" * 60)

                    
            problem = process_event(folder,surrogate_model,results)

            if problem is not None:
                    events_with_issues.append(
                        (folder, problem)
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
