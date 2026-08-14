#main.py
import os
from config import *
from process_event import process_event
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))


def main():

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

        print("\n" + "=" * 60)
        print(f"Processing event {folder}")
        print(f"Event {i}/{total_events}")
        print("=" * 60)

    try:

        problem = process_event(folder)

        if problem is not None:
            events_with_issues.append(
                (folder, problem)
            )

    except ZeroDivisionError as e:

        print(f"Exception while processing {folder}: {e}")

        events_with_issues.append(
            (folder, str(e))
        )

    print("\n" + "=" * 60)
    print("Finished")
    print("=" * 60)

    if events_with_issues:

        print("Events with issues:")

        for folder, reason in events_with_issues:
            print(f"  - {folder}: {reason}")

    else:

        print("All events processed successfully.")

if __name__ == "__main__":
    main()