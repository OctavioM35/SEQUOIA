SEQUOIA is a GW parameter estimation algorithm capable of employing the neuronal networks DANSur, NRHybSur3dq8 and NRSur7dq4 and the approximants IMRPhenomXHM and IMRPhenomXO4a to infer GW signals parameters and plot the obtained posterior distributions alongside those published by the LVK collaboration and uploaded in Zenodo. 
The algorithm processes list of GW data, which are assumed to be named such as:
```text

GWdate              <- Event directory
├── GWdate_l1.h5    <- Event interferometer data
├── GWdate_h1.h5    <- Event interferometer data
├── GWdate_v1.h5    <- Event interferometer data
└── GWdate.h5       <- LVK Zenodo file
```

The priors and parameter estimation configuration values used by SEQUOIA can be customized or they can be automatically loaded from the Zenodo file. This can be configurated in the config.py file.

SEQUOIA follows an installation procedure nearly identical to that of the DANSur surrogate model. Detailed installation instructions are available at the DANSur github repository: https://github.com/osvaldogramaxo/DANSur_22. After following this steps the option install_surrogates must be put on True in config.py to install NRHybSur3dq8 and NRSur7dq4 surrogates models.


