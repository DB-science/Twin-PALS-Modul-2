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

```cpp
// Example code to load gamma streams
loadGammaStream("path/to/511keV_data");
loadGammaStream("path/to/1275keV_data");


![PHS_einstellungen](https://github.com/DB-science/Twin-PALS-Modul-2/assets/102671948/aed6b5af-461e-43e0-ae93-669b3d393922)

