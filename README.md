# SEQUOIA: Sequential Gravitational-Wave Event Parameter Estimation Pipeline

SEQUOIA is a GW parameter estimation algorithm capable of employing the neural networks DANSur, NRHybSur3dq8 and NRSur7dq4 to infer GW signal parameters within DANSur's time of validity (chirp mass > 25 M<sub>⊙</sub>) and compare the resulting posterior distributions with those published by the LVK Collaboration on Zenodo. SEQUOIA may use LVK parameter configuration values if selected.

The algorithm takes as input the paths to directories containing lists of GW data, which are assumed to follow this naming convention:

```text
GWdate              <- Event directory
├── GWdate_l1.h5    <- Event interferometer data
├── GWdate_h1.h5    <- Event interferometer data
├── GWdate_v1.h5    <- Event interferometer data
└── GWdate.h5       <- LVK Zenodo file
```

The priors and parameter estimation configuration values used by SEQUOIA can be customized or automatically loaded from the Zenodo file. This is configured in `config.py`.

## Project structure

```text
SEQUOIA
├── config.py                          # User settings (event, surrogate, paths)
├── main.py                            # Entry point: loops over events
└── sequoia
    │
    ├── utils                          # 1. Setup
    │   ├── generate_dictionary.py     # Saves LVK values used for the inference of the event
    │   └── chirp_mass.py              # Checks precomputed chirp mass
    │
    ├── data                           # 2. Data loading
    │   ├── load_zenodo_directories.py # Locates LVK Zenodo files for the event
    │   └── load_interferometer.py     # Builds detectors (strain, PSD, frequencies)
    │
    ├── priors                         # 3. Priors
    │   ├── load_zenodo_priors.py      # Reads the LVK priors
    │   └── custom_priors.py           # Reads custom priors
    │
    ├── supported_surrogates           # 4. Surrogate selection
    │   ├── surrogates.py              # Registry of available surrogates
    │   ├── install_surrogates.py      # Downloads/installs surrogate data
    │   └── evaluate_selected_surrogate.py  # Calls the chosen surrogate's generator
    │
    ├── waveforms                      # 5. Waveform generation (one per model)
    │   ├── dansur
    │   │   ├── waveform_generator.py  # DANSur waveform generator
    │   │   └── waveform_conversion.py # Parameter conversion for DANSur
    │   ├── nrsur7dq4
    │   │   ├── nrsur_waveform_generator.py  # NRSur7dq4 waveform generator
    │   │   └── nrsur_convert.py       # Parameter conversion for NRSur7dq4
    │   └── nrybsur
    │       ├── nrh_waveform_generator.py    # NRHybSur waveform generator
    │       └── nrh_waveform_conversion.py   # Parameter conversion for NRHybSur
    │
    ├── pipeline                       # 6. Inference
    │   ├── process_event.py           # Runs all steps for one event
    │   └── run_inference_bilby.py     # Bilby + nessai sampling, saves posteriors
    │
    └── plotting                       # 7. Output
        └── plot_manual_corner.py      # Corner plots (surrogate and vs. LVK)
```

## Pipeline

After loading a particular configuration, if selected, SEQUOIA reads the LVK Zenodo file of the event to be inferred and extracts the configuration used by the LVK, in order to apply it in its own parameter estimation.

If a custom configuration is specified in `config.py`, SEQUOIA performs the inference using the provided values.

### Output produced

```text
outdir_NN_GWdate
├── selected_surrogate_corner.png                  # Corner plot of the posteriors
├── selected_surrogate_execution_time.json         # Wall-clock time of the inference
├── selected_surrogate.log                         # Run log
├── selected_surrogate_nessai/                     # nessai sampler output
├── selected_surrogate_result.json                 # Posterior samples
├── selected_surrogate_vs_Zenodo_corner_comparison_GWXXXXXX_XXXXXX.png  # Comparison with LVK
└── GWdate selected_surrogate_reconstructed_waveform.png                # Reconstructed waveform
```

## Setup

SEQUOIA follows an installation procedure nearly identical to that of the DANSur surrogate model. Detailed instructions are available in the [DANSur repository](https://github.com/osvaldogramaxo/DANSur_HM). After following those steps, set `install_surrogates = True` in `config.py` to install the NRHybSur3dq8 and NRSur7dq4 surrogate models.

Run from the parent directory with:

```bash
python3 SEQUOIA/main.py
```

---

# SEQUOIASXS: Sequential Gravitational-Wave Event Parameter Estimation Pipeline for SXS Simulations

> **Under construction**

SEQUOIASXS is a GW parameter estimation algorithm designed with the same philosophy as SEQUOIA, but for SXS simulations. It is nearly identical to SEQUOIA.

## Project structure

```text
SEQUOIASXS
├── config.py                          # User settings (event, surrogate, paths)
├── main.py                            # Entry point: loops over events
└── sequoia
    │
    ├── utils                          # 1. Setup
    │   ├── generate_dictionary.py     # Builds the per-event settings dictionary
    │   └── synthetic_waveform.py      # Converts the SXS waveform into LVK-like data
    │
    ├── data                           # 2. Data loading
    │   ├── load_zenodo_directories.py # Locates LVK Zenodo files for the event
    │   └── load_interferometer.py     # Builds detectors (strain, PSD, frequencies)
    │
    ├── priors                         # 3. Priors
    │   └── custom_priors.py           # Reads custom priors
    │
    ├── supported_surrogates           # 4. Surrogate selection
    │   ├── surrogates.py              # Registry of available surrogates
    │   ├── install_surrogates.py      # Downloads/installs surrogate data
    │   └── evaluate_selected_surrogate.py  # Calls the chosen surrogate's generator
    │
    ├── waveforms                      # 5. Waveform generation (one per model)
    │   ├── dansur
    │   │   ├── waveform_generator.py  # DANSur waveform generator
    │   │   └── waveform_conversion.py # Parameter conversion for DANSur
    │   ├── nrsur7dq4
    │   │   ├── nrsur_waveform_generator.py  # NRSur7dq4 waveform generator
    │   │   └── nrsur_convert.py       # Parameter conversion for NRSur7dq4
    │   └── nrybsur
    │       ├── nrh_waveform_generator.py    # NRHybSur waveform generator
    │       └── nrh_waveform_conversion.py   # Parameter conversion for NRHybSur
    │
    ├── pipeline                       # 6. Inference
    │   ├── process_event.py           # Runs all steps for one event
    │   └── run_inference_bilby.py     # Bilby + nessai sampling, saves posteriors
    │
    └── plotting                       # 7. Output
        └── plot_manual_corner.py      # Corner plots (surrogate and vs. LVK)
```

## Pipeline

After loading a simulation, SEQUOIASXS performs parameter estimation using custom priors and custom configuration values provided through `config.py`.

### Output produced

```text
outdir_NN_SXS_simulation
├── selected_surrogate_corner.png                  # Corner plot of the posteriors
├── selected_surrogate_execution_time.json         # Wall-clock time of the inference
├── selected_surrogate.log                         # Run log
├── selected_surrogate_nessai/                     # nessai sampler output
├── selected_surrogate_result.json                 # Posterior samples
├── Obtained vs inyected parameters for surrogate selected_surrogate.json  # Obtained vs injected values
└── Reconstructed waveform for event SXS:BBH:identifier with selected_surrogate.jpg
```

## Setup

SEQUOIASXS follows the same installation procedure as SEQUOIA (see the [DANSur repository](https://github.com/osvaldogramaxo/DANSur_HM)). After following those steps, set `install_surrogates = True` in `config.py`.

Run from the parent directory with:

```bash
python3 SEQUOIASXS/main.py
```

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.

## Acknowledgements

SEQUOIA includes code adapted from [DANSur](https://github.com/osvaldogramaxo/DANSur_HM),
licensed under the MIT License. Copyright (c) AÑO Autor Original.
