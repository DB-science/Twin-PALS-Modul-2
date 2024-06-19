# Twin-PALS-Modul-2
A combination of a event discriminator and the DLTPulseGenerator to handle the output streams from Modul 1 and create digitized PMT pulses.
This repository contains pre-configured programs for generating positron lifetime spectra. The available components are:

- **all_spectrum**: Complete spectrum
- **ideal_spectrum**: Ideal spectrum
- **Background**: Background noise
- **IRF_ideal**: Ideal Instrument Response Function (IRF)
- **1275_BS-S**: Backscatter of 1275 keV
- **Pile_Up**: Pile-Up events

## Introduction

Positron Annihilation Lifetime Spectroscopy (PALS) is a powerful technique used to study the microstructure of materials. By analyzing the time between the emission of a positron and the detection of the resulting gamma rays, researchers can gain insights into the atomic and electronic structure of the material. This repository provides a suite of tools to simulate various components of a positron lifetime spectrum, allowing users to model and analyze different scenarios.

## Usage

### Preparation

1. **Setup Visual Studio Project**:
    - Integrate all `.cpp` and `.h` files into a Visual Studio project. Ensure all necessary dependencies are included.
    - The main file to modify for your simulation needs is `PALSStudy.cpp`.

2. **Enter Necessary Data**:
    - Within `PALSStudy.cpp`, you'll need to input data such as the gamma energy streams and the parameters for the positron source and lifetime components.

### Data Loading

The first step in your simulation is to load the binary streams for the gamma energies from the Modul 1. These streams represent the detected gamma rays at 511 keV and 1275 keV in your defined digital world.

![Input_Stream_from_modul1](https://github.com/DB-science/Twin-PALS-Modul-2/assets/102671948/852e71a9-c8d4-4d6f-8155-52e3ceac48cb)

### Define PHS limits

Here you define the conversion factor from number of photons on the sensitive area to the pulse height of the PMT pulse in mV. Also you can chose limits for your START and STOP windows in the PHS. Here it's important to define a lower limit for the STOP window to produce minor amount of data. (see picture for example)
![PHS_einstellungen](https://github.com/DB-science/Twin-PALS-Modul-2/assets/102671948/c744483e-0f3a-4803-9d12-cd9a6173210d)

### Parameter Input
After loading the data, you need to set various parameters for the simulation:

- Source Strength: Define the strength of the positron source.
- Sweep: You can chose the time window of your digitizer.
- Pulse Shape: Here you can define the rise time and the pulse width of the resulting PMT pulse.

![define_source_and_pulse_shape](https://github.com/DB-science/Twin-PALS-Modul-2/assets/102671948/0a6def47-ab0f-41ee-96d3-317a10e8abac)

### Define Ground Truth Positron Lifetime

You can define a ground truth (GT) positron lifetime spectrum with up to 5 different lifetime components.
![Define_lifetime_GT](https://github.com/DB-science/Twin-PALS-Modul-2/assets/102671948/5edbcbeb-8210-43d4-8b17-3e55dcf3181b)

### License
This project is licensed under the MIT License.

### Contact
For questions or comments, please contact Your dominik.boras@uni-wuerzburg.de.

# How to cite this Software? 
<p style="font-size:16px;">You must cite the applied version of this software in your study.
<p style="font-size:16px;">You can cite all versions by using the DOI 10.5281/zenodo.11470161. This DOI represents all versions, and will always resolve to the latest one

