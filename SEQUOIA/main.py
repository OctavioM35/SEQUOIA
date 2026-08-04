import os
from config import *
from load_zenodo_directories import load_zenodo
from custom_priors import load_custom_priors
from load_interferometer import load_ifos
from run_inference_bilby import run_inference
from plot_manual_corner import  plot_manual_corner
from chirp_mass import chirp_mass
from load_zenodo_configuration import load_zenodo_configuration
from synthetic_waveform import synthetic_waveform
import bilby
import numpy as np
import sys


def main():
    if generate_synthetic_signal == True:
        print('Generating synthetic waveform...')
        synthetic_waveform(particular_event,output_synthetic)
        print('Synthetic waveform generated succesfuly. Results in: ', output_synthetic)
    else:
        # injection_parameters = {
        #     "chirp_mass": synthetic_chirp_mass,
        #     "mass_ratio": synthetic_mass_ratio,
        #     "spin_1z": synthetic_spin_1,
        #     "spin_2z": synthetic_spin_2,
        #     "luminosity_distance": synthetic_luminosity_distance,
        #     "theta_jn": synthetic_theta_jn,
        #     "ra": synthetic_ra,
        #     "dec": synthetic_declination,
        #     "psi": synthetic_psi,
        #     "phase": synthetic_psi,
        #     "geocent_time": synthetic_gps_time,
        # }
        events_not_supported = []
        folders = [c for c in os.listdir(data_folder) if os.path.isdir(os.path.join(data_folder, c))]
        number_of_events = len(folders)
        for folder in folders:
                # try:
                    number_of_events -=1
                    if run_particular_event is True:
                        if folder != particular_event:
                                continue
                    print(' ')
                    print('=' *50)
                    print(f"Processing event {folder}")
                    if run_particular_event == False:
                        print('Number of events left: ' , number_of_events)
                    event_dir = os.path.join(data_folder, folder)

                    zenodo_file = load_zenodo(data_folder,folder,event_dir)
                    if zenodo_file ==None:
                        print('=' *50)
                        print(f"There is no zenodo data on your folder for this event, {folder}")
                        print('=' *50)
                        events_not_supported.append((folder, 'Zenodo missing'))
                        continue
                    supported_mass = chirp_mass(zenodo_file)
                    if supported_mass == False:
                        print('=' *50)
                        print('DANSur has not been trained for handling this particular event (low or high mass event)')
                        print('=' *50)
                        events_not_supported.append((folder, 'DANSur does not support'))
                        continue
                    outdir = os.path.join(results, f"outdir_NN_{folder}")

                    os.makedirs(outdir, exist_ok=True)

                    if os.path.exists(os.path.join(outdir, 'DANSur_result.json')) and resume == False:
                        print('Event already analyzed')
                        continue

                    if automatic == True:
                        dicc, priors, targ_keys= load_zenodo_configuration(zenodo_file)

                    elif automatic == False:
                        targ_keys ,priors= load_custom_priors()
                    else:
                        print('=' *50)
                        priors('automatic must either be True or False')
                        print('=' *50)
                        sys.exit(1)

                    if targ_keys == None or priors == None:
                            print('=' *50)
                            print('Zenodo priors not valid for this event')
                            print('=' *50)
                            events_not_supported.append((folder, 'Zenodo priors not valid'))
                            continue

                    ifos = load_ifos(event_dir, folder, dicc, zenodo_file,duration)



                    if ifos == None:
                        print('=' *50)
                        print(f"There is no data on your folder for this event, {folder}")
                        print('=' *50)
                        events_not_supported.append((folder, 'Interferometer data missing'))
                        continue
                    


                    run_inference(ifos, dicc,targ_keys, priors, npoints, outdir, resume,approximant, zenodo_file)
                    print('=' *50)
                    print('Plotting results')
                    plot_manual_corner(folder, zenodo_file, event_dir, outdir, approximant)


                    print(f"OK event: {folder}")

                
                # except Exception as e:
                #     print('Excepción', e)
                #     events_not_supported.append((folder,e))
                #     continue
        print('Finished')
        print("Events with issues:  ", events_not_supported)


if __name__ == "__main__":
    main()