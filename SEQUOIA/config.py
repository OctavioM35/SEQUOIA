# config.py

results = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/"
data_folder = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/sintetico"
output_synthetic = '/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/sintetico'
run_particular_event = True             #True for running one event, False for running all events on data_folder directory
particular_event = 'GW150914'
resume = False                          #True for resuming one event if stopped previously, False for running the event from the begining
approximant = False                     #True: IMR estimates waveform , False: DANSur estimates the waveform
automatic = True                        #True for using all Zenodo priors, duration and sampling frecuencies. False for customs
generate_synthetic_signal = False

npoints = 500
duration = 8
sampling_frequency = 1024               #Hz


#Generation of synthetic singal:

synthetic_chirp_mass = 30.69       # masas solares
synthetic_mass_ratio = 1/0.88      # m1/m2

synthetic_gps_time = 1.1263e+9
    # synthetic_geoscent_time = 1126259462.3047757
synthetic_spin_1 = -0.05
synthetic_spin_2 = -0.01

synthetic_luminosity_distance = 473  # Mpc

synthetic_psi = 1.45
synthetic_declination = -1.19         # rad
synthetic_ra = 2.03                   # rad
synthetic_theta_jn = 2.70             # rad

# MIRAR EN DETALLE ESTOS EVENTOS (SON LOS RUIDOSOS)
# ['GW230628_231200', 'GW230819_171910', 'GW231206_233901', 'GW230825_041334', 'GW240109_050431', 'GW231029_111508', 'GW240107_013215', 'GW230814_230901', 'GW230820_212515', 'GW230601_224134', 'GW230803_033412', 'GW231001_140220', 'GW231004_232346', 'GW231113_122623', 'GW230914_111401', 'GW230814_061920']


