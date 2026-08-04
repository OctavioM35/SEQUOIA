import os
import bilby
import h5py
import numpy as np
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
from bilby import likelihood

def load_ifos(event_dir, event, dicc, zenodo_file,duration):

    ifos_list = []
    files = os.listdir(event_dir)
    print('=' *50)
    print('Loading interferometer data...')
    print(' ')
    for archivo in os.listdir(event_dir): 
            nombre, ext = os.path.splitext(archivo)
            if "h1" in nombre:                                                   #Compuebra la que interferómetros están disponibles para el evento concreto
                        h1 = bilby.gw.detector.InterferometerList(["H1"])[0]  # tomar el objeto Interferometer
                        h1_ts = TimeSeries.read(f'{event_dir}/{event}_h1.hdf5' ).resample(int(dicc["sampling-frequency"]))

                        try:
                            with h5py.File(zenodo_file) as file:
                                ap_key = list(file.keys())[0]
                                fs, data = file[f'{ap_key}/psds/H1'][()].T


                            h1_asd = FrequencySeries(data, frequencies=fs)

                            t_gps = h1_ts.times.value[0]+8
                            t_init = t_gps-duration+1
                            t_end = t_init+duration
                            if np.isnan(h1_ts).any():
                                # print('\n\n\n\n', '#'*50,'AHHHHHH', '#'*50, '\n\n\n\n')
                                h1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), h1_ts.duration.value, start_time=h1_ts.times.value[0])
                            else:
                                h1.set_strain_data_from_gwpy_timeseries(time_series=h1_ts.crop(t_init, t_end ))
                                h1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(h1_asd.frequencies, h1_asd.value)
                                ifos_list.append(h1)
                                print('H1 interferometer loaded succesfuly')
                                print(' ')
                        except KeyError:
                              print('No psd found for this interferometer')
                        
                    
            elif "l1" in nombre:
                        l1 = bilby.gw.detector.InterferometerList(["L1"])[0]
                        l1_ts = TimeSeries.read(f'{event_dir}/{event}_l1.hdf5' ).resample(int(dicc["sampling-frequency"]))                        

                        try:

                            with h5py.File(zenodo_file) as file:
                                ap_key = list(file.keys())[0]
                                fs, data = file[f'{ap_key}/psds/L1'][()].T


                            l1_asd = FrequencySeries(data, frequencies=fs)
                            t_gps = l1_ts.times.value[0]+8
                            t_init = t_gps-duration+1
                            t_end = t_init+duration
                            if np.isnan(l1_ts).any():
                                l1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), l1_ts.duration.value, start_time=l1_ts.times.value[0])
                            else:
                                l1.set_strain_data_from_gwpy_timeseries(time_series=l1_ts.crop(t_init, t_end ))

                                l1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(l1_asd.frequencies, l1_asd.value)
                                ifos_list.append(l1)                                
                                print('L1 interferometer loaded succesfuly')
                                print(' ')
                        except KeyError:
                              print('No psd found for this interferometer')
                        
                    
            elif "v1" in nombre.lower():
                            v1 = bilby.gw.detector.InterferometerList(["V1"])[0]
                            v1_ts = TimeSeries.read(f'{event_dir}/{event}_v1.hdf5' ).resample(int(dicc["sampling-frequency"]))

                            try:
                                    with h5py.File(zenodo_file) as file:
                                            ap_key = list(file.keys())[0]
                                            fs, data = file[f'{ap_key}/psds/V1'][()].T
                                            v1_asd = FrequencySeries(data, frequencies=fs)
                                            v1.calibration_model = bilby.gw.detector.calibration.CubicSpline(
                                # envelope_file = file[f"{ap_key}/calibration_envelope/V1"][()],
                                prefix="recalib_V1_",
                                minimum_frequency=int(dicc["minimum-frequency"]["V1"]),
                                maximum_frequency=int(dicc["maximum-frequency"]["V1"]),
                                n_points=int(dicc["spline-calibration-nodes"])
                            )


                                    t_gps = v1_ts.times.value[0]+8
                                    t_init = t_gps-duration+1
                                    t_end = t_init+duration
                                    if np.isnan(v1_ts).any():
                                        v1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), v1_ts.duration.value, start_time=v1_ts.times.value[0])
                                    else:
                                        v1.set_strain_data_from_gwpy_timeseries(time_series=v1_ts.crop(t_init, t_end ))
                                        v1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(v1_asd.frequencies, v1_asd.value)
                                        ifos_list.append(v1)
                                        print('V1 interferometer loaded succesfuly')
                                        print(' ')
                            except KeyError:
                              print('No psd found for this interferometer')

            elif "synthetic" in nombre.lower():

                            # ==========================================================
                            # Detector L1
                            # ==========================================================

                            sin = bilby.gw.detector.InterferometerList(["H1"])[0]

                            # ==========================================================
                            # Leer señal sintética
                            # ==========================================================

                            filename = f'{event_dir}/{event}_synthetic.hdf5'

                            with h5py.File(filename, 'r') as f:
                                strain = f['strain_L1'][:]
                                times_gps = f['times_gps'][:]

                            # ==========================================================
                            # La señal sintética original:
                            # 65536 samples = 8 segundos a 8192 Hz
                            # ==========================================================

                            fs_original = 8192

                            sin_ts = TimeSeries(
                                strain,
                                sample_rate=fs_original,
                                t0=times_gps[0]
                            )

                            # ==========================================================
                            # Información antes del resample
                            # ==========================================================

                            fs = sin_ts.sample_rate.value
                            dt = 1.0 / fs
                            duration_ts = len(sin_ts) / fs

                            print("\nSYNTHETIC SIGNAL")
                            print("---------------------------")
                            print("Samples:", len(sin_ts))
                            print("Sample rate:", fs, "Hz")
                            print("Delta t:", dt, "s")
                            print("Duration:", duration_ts, "s")

                            # ==========================================================
                            # Resample
                            # ==========================================================

                            fs_target = int(dicc["sampling-frequency"])

                            sin_ts = sin_ts.resample(fs_target)

                            # ==========================================================
                            # Información después del resample
                            # ==========================================================

                            fs = sin_ts.sample_rate.value
                            dt = 1.0 / fs
                            duration_ts = len(sin_ts) / fs

                            print("\nAFTER RESAMPLE")
                            print("---------------------------")
                            print("Samples:", len(sin_ts))
                            print("Sample rate:", fs, "Hz")
                            print("Delta t:", dt, "s")
                            print("Duration:", duration_ts, "s")
                            try:
                                    t_gps = sin_ts.times.value[0] + 7.8
                                    t_init = t_gps-duration+1
                                    t_end = t_init+duration
                                    # print(sin_ts,t_init,t_end)
                                    # raise Exception
                                    if np.isnan(sin_ts).any():
                                        sin.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), sin_ts.duration.value, start_time=sin_ts.times.value[0])
                                    else:
                                        sin.set_strain_data_from_gwpy_timeseries(time_series=sin_ts.crop(t_init, t_end ))
                                        print(sin_ts)
                                        print(len(sin_ts))
                                        ifos_list.append(sin)
                                        # raise Exception                                
                                        print('L1 synthetic interferometer loaded succesfuly')
                                        print(' ')
                                        # raise Exception
                            except KeyError:
                                        print('No psd found for this interferometer')
            
    if len(ifos_list) ==0:
          return None
    print(ifos_list)
    print('Interferometer data loaded succesfully')
    print('=' *50)
    return bilby.gw.detector.InterferometerList(ifos_list)