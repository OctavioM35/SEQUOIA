An algorithm which uses the neuronal network DANSur for estimate GW signals parameters and plots the distributions obtained alongside with the results from Zenodo. 
It works for a list of events, which are assumed to be named such as:

```text

GWdate              <- Event directory
├── GWdate_l1.h5   
├── GWdate_h1.h5    <- Event interferometer data
├── GWdate_v1.h5
└── GWdate.h5       <- Zenodo file
```

These events were obtained used the python gwosc package and will not work for the events from GWTC webpage (though the script can easily be modified to make it work) Right now it usses the priors from the zenodo file but this can be configurated at will on the zenodo priors file. The duration, sampling frecuency and route to the directories with data are configurated in the configuration section. The instalation of DANSur is needed to make it work, instructrions for this can be found in https://github.com/osvaldogramaxo/DANSur_22/

Total pacakages needed: bilby h5py matplotlib numpy os gwpy.timeseries gwpy.frequencyseries bilby.gw.likelihood plus all those needed to make dansur work (see previos link for reference) Keep in mind that DANSur does not currently work for low masses events (chirp mass <20 solar masses).
