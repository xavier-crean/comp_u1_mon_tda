# Topological Data Analysis of Monopoles in U(1) Lattice Gauge Theory - Monte Carlo and Analysis Code Release

## Contents

- [Topological Data Analysis of Monopoles in U(1) Lattice Gauge Theory - Monte Carlo and Analysis Code Release](#topological-data-analysis-of-monopoles-in-u1-lattice-gauge-theory---monte-carlo-and-analysis-code-release)
  - [Contents](#contents)
  - [Overview](#overview)
  - [Code to generate lattice configuration data (Step 1)](#code-to-generate-lattice-configuration-data-step-1)
  - [Code to process and analyse configuration data (Step 2 \& Step 3)](#code-to-process-and-analyse-configuration-data-step-2--step-3)
    - [Setup](#setup)
      - [Requirements](#requirements)
      - [Installation](#installation)
      - [Data Release](#data-release)
    - [Binder](#binder)
    - [Processing the raw lattice configuration data (Step 2)](#processing-the-raw-lattice-configuration-data-step-2)
    - [Performing the analysis (Step 3)](#performing-the-analysis-step-3)
      - [Cached data files](#cached-data-files)
  - [Project organisation](#project-organisation)

## Overview

This repository accompanies the paper [Topological Data Analysis of Monopoles in U(1) Lattice Gauge Theory][paper]. It contains code for

1. [generating U(1) lattice gauge field configurations][step_1] via Monte Carlo simulation, 
2. [computing observables][step_2] from lattice configurations, 
3. [producing plots and results][step_3] included in the [paper][paper].

A [directory tree](#project-organisation) is included at the bottom of this file.

Alongside this repository, there exists [an accompanying data release][data] where [Step 1][step_1] and [Step 2][step_2] in the pipeline have been pre-computed. Therefore, to reproduce analysis from the paper using the data release, the user is referred to [Step 3][step_3] ([setup instructions](#setup) are included below).

Using Binder, a Docker container has been constructed for reproduction of results from the paper (without the need to install Python or dependencies) and can be accessed online ([see below](#binder)).

## Code to generate lattice configuration data (Step 1)

The code to generate U(1) lattice gauge field configurations, which is largely written in C, forms part of a wider package for simulating lattice gauge field theories. It must be configured, built and run using a Linux distribution. The software is stored in `src/data/` and contains its own

* [README][mc_README] (software details),
* [INSTALL][mc_INSTALL] file (installation instructions)
* and [AUTHORS][mc_AUTHORS] file.

Using this code, one may reproduce the ensemble averages quoted in the [paper][paper]. Note that for bit-for-bit reproducibility one must set a fixed random seed and run using the same distribution each time.

## Code to process and analyse configuration data (Step 2 & Step 3)

The code to perform analysis of lattice gauge field configurations is written in Python. See Section 2.2 and Section 3 in the [paper][paper]. It consists of two steps:

* [Step 2][step_2]: Raw lattice configuration data is processed by the scripts in directory `src/observables/`. Collectively, this may take several hours to compute. For reproducibility, our processed data is included in the [data release][data] and, once extracted, is stored in directory `data/observables/`.
* [Step 3][step_3]: Plots and results from the [paper][paper] may be reproduced by running the script `src/analysis/analysis.py` which takes about an hour to run. This will generate [cached data](#cached-data-files), stored in `cached_data/`, for quicker re-runs. 
  
  **N.B.** Pre-computed cached data has also been provided in the [data release][data]; with cached data, `src/analysis/analysis.py` takes about one minute to run.

In this section of the pipeline, data files are stored in the [HDF5 format][hdf5].

---

### Setup

#### Requirements

The code has been tested with Python 3.10. 

Dependencies are documented in `environment.yml` and are most easily managed via a Conda environment. 

The commands below refer to a Linux or compatible environment and need to be run from the root directory of the repository.

#### Installation

* Download the repository
* Create a new Conda environment with the necessary dependencies using

      conda env create -f environment.yml
    (Alternatively, create a Python environment with the listed packages.)

#### Data Release

* Download the `data.tar.gz` from [the accompanying data release][data] and extract into the root directory of the repository using

        tar -xf data.tar.gz
* To use **pre-computed** [cached data files](#cached-data-files), for fast reproduction of the results in the paper, download `cached_data.tar.gz` from the [data release][data] and extract into the root directory of the repository using

        tar -xf cached_data.tar.gz

---

### Binder

Using [mybinder.org][binder], a Docker container with [necessary dependencies](#requirements) has been constructed for reproduction of the results in the paper. Note that the [data release][data] must be uploaded and extracted into the root directory (see [data release setup](#data-release)).

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/xavier-crean/comp_u1_mon_tda.git/HEAD)

---

### Processing the raw lattice configuration data (Step 2)

The directory `src/observables/` contains two python scripts and a python module:

* `src/observables/action.py` for computing the total action of a lattice configuration; it takes four command line arguments: lattice size L, $\beta$ value (to 4d.p.), number of samples and number of sweeps between measurement. E.g.,

      python src/observables/action.py 6 0.9000 200 350000
    The output is saved into directory `data/observables/action/L.L.L.L/` in `.h5` format.
* `src/observables/persistence_diagram.py` for computing the zeroth and first order homology of monopole current networks of a lattice configuration; it takes five command line arguments: lattice size L, $\beta$ (to 4d.p.), number of samples, number of sweeps between measurement and number of parallel computations. E.g.,

      python src/observables/persistence_diagram.py 6 0.9000 200 350000 2
    The output is saved into directory `data/observables/pd/L.L.L.L/` in `.h5` format.
* `src/observables/configurations.py` a module used to parse lattice configurations from raw configuration data files.

---

### Performing the analysis (Step 3)

The python script `src/analysis/analysis.py` which

* makes directories `cached_data/`, `reports/` and `reports/figures/` (if they don't exist already),
* computes the average action observable and the Betti number observables, then saves figures in `.pdf` format into `reports/figures/`,
* and outputs results from a finite-size scaling analysis into `reports/`
  * in `.csv` format (with header)
  * and in `.tex` format (for a LaTeX table).

The script is run using

      python src/analysis/analysis.py

#### Cached data files

Intermediary files, for caching and running the code quickly second time round, are stored in the directory `cached_data/`. To use the files in `cached_data.tar.gz` (from the [data release][data]), it is important to ensure the correct random seeds are used. With the cached data files, `src/analysis/analysis.py` takes less than a minute to run; without them, it takes about an hour.

**N.B.** if `cached_data.tar.gz` is extracted into the root, then `src/analysis/analysis.py` will read from this directory. If you are trying to reproduce the analysis using `data.tar.gz`, you will need to delete or rename the directory `cached_data/`.

---

## Project organisation

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── observables    <- Action and persistence diagrams saved in .h5 format.
    |   └── configurations <- Lattice gauge field configuration data stored in conf.dat in IDLG-like format.
    │
    ├── cached_data        <- Intermediary cached data files in .h5 format for faster reproduction of results.
    ├── reports            <- Generated tables used in the paper as .csv and .tex.
    │   └── figures        <- Generated figures used in the paper as .pdf.
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

[data]: https://doi.org/10.5281/zenodo.10806046
[paper]: https://arxiv.org/abs/2403.07739
[mc_README]: src/data/README
[mc_AUTHORS]: src/data/AUTHORS
[mc_install]: src/data/INSTALL
[binder]: https://mybinder.org/
[hdf5]: https://www.hdfgroup.org/solutions/hdf5
[step_1]: #code-to-generate-lattice-configuration-data-step-1
[step_2]: #processing-the-raw-lattice-configuration-data-step-2
[step_3]: #performing-the-analysis-step-3