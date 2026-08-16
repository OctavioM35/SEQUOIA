import os
import bilby
import h5py
import numpy as np
from gwpy.frequencyseries import FrequencySeries
from gwpy.timeseries import TimeSeries
from scipy.interpolate import interp1d


def load_ifos(event_dir, event, dicc, zenodo_file):
    duration=dicc["duration"]
    sampling_frequency =dicc["sampling-frequency"]
    ifos_list = []
    print('=' *50)
    print('Loading interferometer data...')
    print(' ')
    for folder in os.listdir(event_dir): 
            name, ext = os.path.splitext(folder)

            if "h1" in name:                                                   
                        h1 = bilby.gw.detector.InterferometerList(["H1"])[0]  
                        h1_ts = TimeSeries.read(f'{event_dir}/{event}_h1.hdf5' ).resample(int(dicc["sampling-frequency"]))
                        h1.minimum = dicc["frequencies"]["h1"]["minimum"]
                        h1.maximun = dicc["frequencies"]["h1"]["maximum"]
                        try:

                            # ---------------------------------------------------------
                            # 7.1 Load Zenodo PSD
                            # ---------------------------------------------------------  
                            
                            with h5py.File(zenodo_file) as file:
                                ap_key = list(file.keys())[0]
                                fs, data = file[f'{ap_key}/psds/H1'][()].T
                            h1_asd = FrequencySeries(data, frequencies=fs)

                            t_gps = dicc["gps_time"]
                            t_init = t_gps-duration+1
                            t_end = t_init+duration

                            if np.isnan(h1_ts).any():
                                print('------' *10)
                                print('Warning: missing data on interferometer h1 found. Skipping this interferometer...')
                                h1.set_strain_data_from_zero_noise(sampling_frequency, h1_ts.duration.value, start_time=h1_ts.times.value[0])
                            else:
                                h1.set_strain_data_from_gwpy_timeseries(time_series=h1_ts.crop(t_init, t_end ))
                                h1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(h1_asd.frequencies, h1_asd.value)
                                ifos_list.append(h1)
                        
                        except KeyError:
                              print('Warning: No psd found for interferometer H1. Skipping interferometer...')
                        
            elif "l1" in name:
                        l1 = bilby.gw.detector.InterferometerList(["L1"])[0]
                        l1_ts = TimeSeries.read(f'{event_dir}/{event}_l1.hdf5' ).resample(int(dicc["sampling-frequency"]))                        
                        l1.minimum = dicc["frequencies"]["l1"]["minimum"]
                        l1.maximun = dicc["frequencies"]["l1"]["maximum"]
                        try:

                            with h5py.File(zenodo_file) as file:
                                ap_key = list(file.keys())[0]
                                fs, data = file[f'{ap_key}/psds/L1'][()].T


                            l1_asd = FrequencySeries(data, frequencies=fs)
                            t_gps = dicc["gps_time"]
                            t_init = t_gps-duration+1
                            t_end = t_init+duration
                            if np.isnan(l1_ts).any():
                                print('------' *10)
                                print('Warning: missing data on interferometer l1. Skipping this interferometer... ')

                                l1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), l1_ts.duration.value, start_time=l1_ts.times.value[0])
                            else:
                                l1.set_strain_data_from_gwpy_timeseries(time_series=l1_ts.crop(t_init, t_end ))

                                l1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(l1_asd.frequencies, l1_asd.value)
                                ifos_list.append(l1)                                
                                print('L1 interferometer loaded succesfuly')
                                print(' ')
                        except KeyError:
                              print('Warning: No psd found for interferometer L1. Skipping interferometer...')
                        

                                            
            elif "v1" in name.lower():
                        v1 = bilby.gw.detector.InterferometerList(["V1"])[0]
                        v1_ts = TimeSeries.read(f'{event_dir}/{event}_v1.hdf5' ).resample(int(dicc["sampling-frequency"]))
                        v1.minimum = dicc["frequencies"]["v1"]["minimum"]
                        v1.maximun = dicc["frequencies"]["v1"]["maximum"]
                        try:
                                    with h5py.File(zenodo_file) as file:
                                            ap_key = list(file.keys())[0]
                                            fs, data = file[f'{ap_key}/psds/V1'][()].T
                                            v1_asd = FrequencySeries(data, frequencies=fs)
                                #             v1.calibration_model = bilby.gw.detector.calibration.CubicSpline(
                                # # envelope_file = file[f"{ap_key}/calibration_envelope/V1"][()],
                                # prefix="recalib_V1_",
                                # minimum_frequency=int(dicc["minimum-frequency"]["V1"]),
                                # maximum_frequency=int(dicc["maximum-frequency"]["V1"]),
                                # n_points=int(dicc["spline-calibration-nodes"])
                            # )

                                    t_gps = dicc["gps_time"]
                                    t_init = t_gps-duration+1
                                    t_end = t_init+duration
                                    if np.isnan(v1_ts).any():
                                        print('------' *10)
                                        print('Warning: missing data on interferometer v1. Skipping this interferometer...')
                                        v1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), v1_ts.duration.value, start_time=v1_ts.times.value[0])
                                    else:
                                        v1.set_strain_data_from_gwpy_timeseries(time_series=v1_ts.crop(t_init, t_end ))
                                        v1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(v1_asd.frequencies, v1_asd.value)
                                        ifos_list.append(v1)
                                        print('V1 interferometer loaded succesfuly')
                                        print(' ')
                        except KeyError:
                              print('Warning: No psd found for interferometer V1. Skipping interferometer...')
        
            
    if len(ifos_list) ==0:
          return None

    print('Interferometer data loaded succesfully')
    print('=' *50)
    return bilby.gw.detector.InterferometerList(ifos_list)