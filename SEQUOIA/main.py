import os
from config import *
from load_zenodo_directories import load_zenodo
from zenodo_priors import load_zenodo_priors
from load_interferometer import load_ifos
from run_inference_bilby import run_inference
from plot_manual_corner import  plot_manual_corner
import bilby



def main():
    events_not_supported = []
    folders = [c for c in os.listdir(data_folder) if os.path.isdir(os.path.join(data_folder, c))]
    number_of_events = len(folders)
    for folder in folders:
                number_of_events -=1
                if run_particular_event is True:
                    if folder != particular_event:
                            continue
                print(' ')
                print('=' *50)
                print(f"Processing event {folder}")
                print('Number of events left: ' , number_of_events)
                event_dir = os.path.join(data_folder, folder)
                outdir = os.path.join(results, f"outdir_NN_{folder}")

                zenodo_file = load_zenodo(data_folder,folder,event_dir)
                if zenodo_file ==None:
                    print(f"There is no zenodo data on your folder for this event, {folder}")
                    events_not_supported.append((folder, 'Zenodo missing'))
                    continue


                targ_keys ,priors= load_zenodo_priors(zenodo_file)
                if targ_keys == None and priors == None:
                      print('Zenodo priors not valid for this event')
                      events_not_supported.append((folder, 'Zenodo priors not valid'))
                      continue

                ifos = load_ifos(event_dir, folder, sampling_frequency, zenodo_file,duration)
                if ifos == None:
                     print(f"There is no data on your folder for this event, {folder}")
                     events_not_supported.append((folder, 'Interferometer data missing'))
                     continue
                
                run_inference(ifos, duration, sampling_frequency,targ_keys, priors, npoints, outdir, resume)
                print('=' *50)
                print('Plotting results')
                plot_manual_corner(folder, zenodo_file, event_dir, outdir)


                print(f"OK event: {folder}")

    print('Finished')
    print("Events with issues:  ", events_not_supported)



if __name__ == "__main__":
    main()