import os
import bilby
import h5py
import numpy as np
from gwpy.frequencyseries import FrequencySeries
from gwpy.timeseries import TimeSeries
from scipy.interpolate import interp1d

# ---------------------------------------------------------
# 7. Load interferometers
# ---------------------------------------------------------


def load_ifos(event_dir, event, dicc, zenodo_file):
    duration=dicc["duration"]
    ifos_list = []
    print('=' *50)
    print('Loading interferometer data...')
    print(' ')
    for folder in os.listdir(event_dir): 
            name, ext = os.path.splitext(folder)

            if "h1" in name.lower() :                                                   


                h1 = bilby.gw.detector.InterferometerList(["H1"])[0]

                # ---------------------------------------------------------
                # 1. Load strain data
                # ---------------------------------------------------------
                h1_ts = TimeSeries.read(
                    f'{event_dir}/{event}_h1.hdf5'
                ).resample(int(dicc["sampling-frequency"]))

                # ---------------------------------------------------------
                # 2. Frequency range
                # ---------------------------------------------------------
                if dicc["frequencies"]["h1"]["minimum"] != 0:
                    h1.minimum_frequency = dicc["frequencies"]["h1"]["minimum"]

                if dicc["frequencies"]["h1"]["maximum"] != 0:
                    h1.maximum_frequency = dicc["frequencies"]["h1"]["maximum"]

                try:

                    # ---------------------------------------------------------
                    # 3. Load Zenodo PSD
                    # ---------------------------------------------------------
                    with h5py.File(zenodo_file) as file:
                        ap_key = list(file.keys())[0]

                        fs, data = file[
                            f'{ap_key}/psds/H1'
                        ][()].T

                    h1_psd = FrequencySeries(
                        data,
                        frequencies=fs
                    )

                    # ---------------------------------------------------------
                    # 4. Start/end time
                    # ---------------------------------------------------------
                    if dicc["start_time"] == 0:

                        try:
                            dicc["start_time"] = (
                                dicc["trigger_time"]
                                - duration
                                + dicc["post_trigger"]
                            )

                        except Exception:
                            dicc["gps_time"] = (
                                h1_ts.times.value[0] + 8
                            )

                            dicc["start_time"] = (
                                dicc["gps_time"]
                                - duration
                                + 2
                            )

                    t_init = dicc["start_time"]
                    t_end = t_init + duration

                    # ---------------------------------------------------------
                    # 5. Check NaNs
                    # ---------------------------------------------------------
                    if np.isnan(h1_ts).any():

                        print('------' * 10)
                        print(
                            'Warning: missing data on interferometer '
                            'H1. Skipping this interferometer...'
                        )

                        h1.set_strain_data_from_zero_noise(
                            int(dicc["sampling-frequency"]),
                            h1_ts.duration.value,
                            start_time=h1_ts.times.value[0]
                        )

                    else:

                        # -----------------------------------------------------
                        # 6. Load cropped strain into Bilby
                        # -----------------------------------------------------
                        h1.set_strain_data_from_gwpy_timeseries(
                            time_series=h1_ts.crop(
                                t_init,
                                t_end
                            )
                        )

                        # -----------------------------------------------------
                        # 7. Tukey roll-off
                        # -----------------------------------------------------
                        h1.strain_data.roll_off = dicc["roll_off"]

                        # -----------------------------------------------------
                        # 8. PSD
                        # -----------------------------------------------------
                        h1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
                                h1_psd.frequencies,
                                h1_psd.value
                            )
                        

                        # -----------------------------------------------------
                        # 10. Add H1
                        # -----------------------------------------------------
                        ifos_list.append(h1)

                        print('H1 interferometer loaded successfully')
                        print(' ')
                        
                except KeyError:

                    print(
                        'Warning: No PSD found for interferometer H1. '
                        'Skipping interferometer...'
                    )
                                                
            elif "l1" in name.lower() :
                        l1 = bilby.gw.detector.InterferometerList(["L1"])[0]
                        l1_ts = TimeSeries.read(f'{event_dir}/{event}_l1.hdf5' ).resample(int(dicc["sampling-frequency"]))                        

                        if dicc["frequencies"]["l1"]["minimum"] !=0:
                            l1.minimum_frequency = dicc["frequencies"]["l1"]["minimum"]
                        if dicc["frequencies"]["l1"]["maximum"] !=0:
                            l1.maximum_frequency = dicc["frequencies"]["l1"]["maximum"] 

                        try:

                            # ---------------------------------------------------------
                            # 7.2 Load Zenodo PSD and set strain for L1 interferometer
                            # ---------------------------------------------------------  

                            with h5py.File(zenodo_file) as file:
                                ap_key = list(file.keys())[0]
                                fs, data = file[f'{ap_key}/psds/L1'][()].T


                            l1_psd = FrequencySeries(data, frequencies=fs)
                            if dicc["start_time"] ==0:
                                    dicc["gps_time"] = l1_ts.times.value[0]+8
                                    dicc["start_time"] = dicc["gps_time"]-duration+2
                            
                            t_init = dicc["start_time"]
                            t_end = t_init+duration

                            if np.isnan(l1_ts).any():
                                print('------' *10)
                                print('Warning: missing data on interferometer L1. Skipping this interferometer... ')

                                l1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), l1_ts.duration.value, start_time=l1_ts.times.value[0])
                            else:
                                l1.set_strain_data_from_gwpy_timeseries(time_series=l1_ts.crop(t_init, t_end ))
                                l1.strain_data.roll_off = dicc["roll_off"]

                                l1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(l1_psd.frequencies, l1_psd.value)
                                ifos_list.append(l1)                                
                                print('L1 interferometer loaded succesfuly')
                                print(' ')
                        except KeyError:
                              print('Warning: No psd found for interferometer L1. Skipping interferometer...')
                        

                                            
            elif "v1" in name.lower():
                        v1 = bilby.gw.detector.InterferometerList(["V1"])[0]
                        v1_ts = TimeSeries.read(f'{event_dir}/{event}_v1.hdf5' ).resample(int(dicc["sampling-frequency"]))

                        if dicc["frequencies"]["v1"]["minimum"] !=0:
                            v1.minimum_frequency = dicc["frequencies"]["v1"]["minimum"]
                        if dicc["frequencies"]["v1"]["maximum"] !=0:
                            v1.maximum_frequency = dicc["frequencies"]["v1"]["maximum"] 
                        try:
                                    
                            # ---------------------------------------------------------
                            # 7.3 Load Zenodo PSD and set strain for V1 interferometer
                            # ---------------------------------------------------------  
                            
                                    with h5py.File(zenodo_file) as file:
                                            ap_key = list(file.keys())[0]
                                            fs, data = file[f'{ap_key}/psds/V1'][()].T
                                    v1_psd = FrequencySeries(data, frequencies=fs)
                                #             v1.calibration_model = bilby.gw.detector.calibration.CubicSpline(
                                # # envelope_file = file[f"{ap_key}/calibration_envelope/V1"][()],
                                # prefix="recalib_V1_",
                                # minimum_frequency=int(dicc["minimum-frequency"]["V1"]),
                                # maximum_frequency=int(dicc["maximum-frequency"]["V1"]),
                                # n_points=int(dicc["spline-calibration-nodes"])
                            # )

                                    if dicc["start_time"] ==0:
                                            dicc["gps_time"] = v1_ts.times.value[0]+8
                                            dicc["start_time"] = dicc["gps_time"]-duration+2
                                    
                                    t_init = dicc["start_time"]
                                    t_end = t_init+duration

                                    if np.isnan(v1_ts).any():
                                        print('------' *10)
                                        print('Warning: missing data on interferometer V1. Skipping this interferometer...')
                                        v1.set_strain_data_from_zero_noise(int(dicc["sampling-frequency"]), v1_ts.duration.value, start_time=v1_ts.times.value[0])
                                    else:
                                        v1.set_strain_data_from_gwpy_timeseries(time_series=v1_ts.crop(t_init, t_end ))
                                        v1.strain_data.roll_off = dicc["roll_off"]

                                        v1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(v1_psd.frequencies, v1_psd.value)
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