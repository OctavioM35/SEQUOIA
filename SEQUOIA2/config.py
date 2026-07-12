# config.py
duration = 8 
sampling_frequency = 1024               #Hz
npoints = 500
run_particular_event = True             #True for running one event, False for running all events on data_folder directory
particular_event = 'GW230824_033047'
resume = False                          #True for resuming one event if stopped previously, False for running the event from the begining
approximant = False                     #True: IMR estimates waveform , False: DANSur estimates the waveform


results = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master"
data_folder = "/home/octaio-m/PycharmProjects/PythonProject/TFG/DANSur_22-master/gwtc_4"


