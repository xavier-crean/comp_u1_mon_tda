# Topological Data Analysis of Monopoles in U(1) Lattice Gauge Theory

## Overview
This repository contains code for
1. generating U(1) lattice gauge field configurations via Monte Carlo simulation, 
2. computing observables from lattice configurations, 
3. producing plots and results included in the paper [Topological Data Analysis of Monopoles in U(1) Lattice Gauge Theory][paper].

A directory tree is included at the bottom of this file.

## Code to generate lattice configuration data
The code to generate U(1) lattice gauge field configurations, which is largely written in C, forms part of a wider package for simulating lattice gauge field theories. It must be configured, built and run using a Linux distribution. The software is stored in `/src/data/` and contains its own 
* [README.md][mc_README] (software details), 
* [INSTALL][mc_INSTALL] file (installation instructions) 
* and [AUTHORS][mc_AUTHORS] file.

Using this code, one may reproduce the ensemble averages quoted in the [paper][paper]. Note that for bit-for-bit reproducibility one must set a fixed random seed and run using the same distribution each time.

## Code to perform analysis of data
The code to perform analysis of lattice gauge field configurations is written in Python. See Section 2.2 and Section 3 in the [paper][paper].

Raw lattice configuration data is processed by the scripts in directory `/src/observables/`. Collectively, this may take several hours to compute. For reproducibility, our processed data is included in the [data release][data] to be stored in directory `/data/observables/`.

With intermediate saved data from the [data release][data], one may quickly reproduce the results from the [paper][paper] by running the script `/src/analysis/analysis.py`. 

### Setup

#### Requirements

The code has been tested with Python 3.10. 

Dependencies are documented in `environment.yml` and are most easily managed via a Conda environment. 

#### Installation
* Download the repository
* Create a new Conda environment with the necessary dependencies using

      conda env create -f environment.yml
    (Alternatively, create a Python environment with the listed packages.)

#### Data

* Download the data from [the accompanying data release][data] into the directory `/data/`.
* To use intermediate saved files, extract the contents of `intermediate.zip` into the directory `/data/intermediate/`.
---

### Processing the raw lattice configuration data
The directory `/src/observables/` contains two python scripts and a python module:
* `action.py` for computing the total action of a lattice configuration; it takes four command line arguments: lattice size L, $\beta$ value (to 4d.p.), number of samples and number of sweeps between measurement. E.g.,

      python action.py 6 0.9000 200 350000
    The output is saved into directory `/data/observables/action/L.L.L.L/` in `.h5` format.
* `persistence_diagram.py` for computing the zeroth and first order homology of monopole current networks of a lattice configuration; it takes five command line arguments: lattice size L, $\beta$ (to 4d.p.), number of samples, number of sweeps between measurement and number of parallel computations. E.g.,

      python persistence_diagram.py 6 0.9000 200 350000 2
    The output is saved into directory `/data/observables/pd/L.L.L.L/` in `.h5` format.
* `configurations.py` a module used to parse lattice configurations from raw configuration data files.

### Performing the analysis
The directory `/src/analysis/` contains the python script `analysis.py` which
* makes directories `/reports/` and `/reports/figures/` (if they don't exist already),
* computes the average action observable and the Betti number observables, then saves figures in `.png` format into `/reports/figures/`,
* and outputs results from a finite-size scaling analysis into `/reports/`
    * in `.csv` format (with and without header) 
    * and in `.tex` format (for a LaTeX table).

Intermediate saved files, for checkpointing and running the code quickly second time round, are stored in the directory `/data/intermediate/`. To use the files in `intermediate.zip` (from the [data release][data]), it is important to ensure the correct random seeds are used. With the intermediate files, `analysis.py` takes a few minutes to run; without them, it takes a few hours.

The script is run using

      python analysis.py

---
## Project Organisation

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── intermediate   <- Intermediate saved files in .h5 format.
    │   ├── observables    <- Action and persistence diagrams saved in .h5 format.
    |   └── configurations <- Lattice gauge field configuration data stored in conf.dat in IDLG-like format.
    │
    ├── reports            <- Generated tables used in the paper as .csv and .tex.
    │   └── figures        <- Generated figures used in the paper as .png.
    │
    ├── environment.yml    <- The requirements file for reproducing the Conda environment.
    │                         
    └── src                <- Source code for this project.
        │
        ├── data           <- Source code written in C for generating U(1) lattice gauge field 
        |                     configurations via Monte Carlo Markov Chain simulation.
        │
        ├── observables    <- Python scripts to compute action and Betti numbers from configurations
        |   |                 (via a trivial cubical filtration)
        │   ├── action.py
        │   ├── configuration.py
        │   └── persistence_diagram.py   
        │
        └── analysis       <- Python script to plot observables against beta, perform finite-size 
            |                 scaling analysis and generate tables used in the paper.                    
            └── analysis.py
---

[data]: https://doi.org/10.5281/zenodo.7060073
[paper]: https://arxiv.org/abs/2207.13392
[mc_README]: /src/data/README.md,
[mc_AUTHORS]: /src/data/AUTHORS
[mc_install]: /src/data/INSTALL