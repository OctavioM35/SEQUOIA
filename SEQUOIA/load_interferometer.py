import os
import bilby
import h5py
import numpy as np
from gwpy.timeseries import TimeSeries
from gwpy.frequencyseries import FrequencySeries
from bilby import likelihood

def load_ifos(event_dir, event, sampling_frequency, zenodo_file,duration):

    ifos_list = []
    files = os.listdir(event_dir)
    print(event_dir)
    for archivo in os.listdir(event_dir): 
            nombre, ext = os.path.splitext(archivo)
            if "h1" in nombre:                                                   #Compuebra la que interferómetros están disponibles para el evento concreto
                        h1 = bilby.gw.detector.InterferometerList(["H1"])[0]  # tomar el objeto Interferometer
                        h1_ts = TimeSeries.read(f'{event_dir}/{event}_h1.hdf5' ).resample(sampling_frequency)
                        with h5py.File(zenodo_file) as file:
                            ap_key = list(file.keys())[0]
                            fs, data = file[f'{ap_key}/psds/H1'][()].T
                        h1_asd = FrequencySeries(data, frequencies=fs)

                        t_gps = h1_ts.times.value[0]+8
                        t_init = t_gps-duration+1
                        t_end = t_init+duration
                        if np.isnan(h1_ts).any():
                            #print('\n\n\n\n', '#'*50,'AHHHHHH', '#'*50, '\n\n\n\n')
                            h1.set_strain_data_from_zero_noise(sampling_frequency, h1_ts.duration.value, start_time=h1_ts.times.value[0])
                        else:
                            h1.set_strain_data_from_gwpy_timeseries(time_series=h1_ts.crop(t_init, t_end ))
                            h1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(h1_asd.frequencies, h1_asd.value)
                            ifos_list.append(h1)
            elif "l1" in nombre:
                        l1 = bilby.gw.detector.InterferometerList(["L1"])[0]
                        l1_ts = TimeSeries.read(f'{event_dir}/{event}_l1.hdf5' ).resample(sampling_frequency)

                        with h5py.File(zenodo_file) as file:
                            ap_key = list(file.keys())[0]
                            fs, data = file[f'{ap_key}/psds/L1'][()].T
                        l1_asd = FrequencySeries(data, frequencies=fs)
                        t_gps = l1_ts.times.value[0]+8
                        t_init = t_gps-duration+1
                        t_end = t_init+duration
                        if np.isnan(l1_ts).any():
                            #print('\n\n\n\n', '#'*50,'AHHHHHH', '#'*50, '\n\n\n\n')
                            l1.set_strain_data_from_zero_noise(sampling_frequency, l1_ts.duration.value, start_time=l1_ts.times.value[0])
                        else:
                            l1.set_strain_data_from_gwpy_timeseries(time_series=l1_ts.crop(t_init, t_end ))
                            l1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(l1_asd.frequencies, l1_asd.value)
                            ifos_list.append(l1)
                        
                    
            elif "v1" in nombre.lower():
                            v1 = bilby.gw.detector.InterferometerList(["V1"])[0]
                            v1_ts = TimeSeries.read(f'{event_dir}/{event}_v1.hdf5' ).resample(sampling_frequency)
                            t_gps = v1_ts.times.value[0]+8
                            t_init = t_gps-duration+1
                            t_end = t_init+duration
                            v1_asd = v1_ts.psd(fftlength=2)
                            if np.isnan(v1_ts).any():
                                v1.set_strain_data_from_zero_noise(sampling_frequency, v1_ts.duration.value, start_time=v1_ts.times.value[0])
                            else:
                                v1.set_strain_data_from_gwpy_timeseries(time_series=v1_ts.crop(t_init, t_end ))
                                v1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(v1_asd.frequencies, v1_asd.value)
                                ifos_list.append(v1)
            
    if len(ifos_list) ==0:
          return None
    return bilby.gw.detector.InterferometerList(ifos_list)