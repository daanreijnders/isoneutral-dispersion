# Simulating Lagrangian Subgrid-Scale Dispersion on Neutral Surfaces in the Ocean
This repository holds the data generation and analysis scripts for the manuscript 'Simulating Lagrangian Subgrid-Scale Dispersion on Neutral Surfaces in the Ocean'. 
https://github.com/daanreijnders/isoneutral-dispersion
A snapshot of this repository is stored at https://doi.org/10.24416/UU01-RXA2PB

The corresponding publication is available at.
Reijnders, D., Deleersnijder, E., & Sebille, E. (2022). Simulating Lagrangian Subgrid‐Scale Dispersion on Neutral Surfaces in the Ocean. Journal of Advances in Modeling Earth Systems, 14(2). https://doi.org/10.1029/2021MS002850


## `ACC_mitgcm_config`
Contains MITgcm configuration files for the coarse and fine model, including their spinups. This includes the compilation code and namelists.

## `ACC_particle_experiments`
Contains Python scripts for the Lagrangian simulations carried out with [Parcels](https://github.com/OceanParcels/parcels) (experiments with a tracer patch, and with particles initialized on a lattice), and analysis scripts pertaining to these simulations and their output.

## `data_processing`
Contains scripts and notebooks for processing the data. This includes creating coarsened fields, computing discrete derivatives, and computing the quantities necessary for the LS parameterization.

## `idealized_experiments`
Experiments using the idealized set-up. 

## `kernels`
Functions used by [Parcels](https://github.com/OceanParcels/parcels) to drive the Markov models.

## `misc`
Contains a PDF copy of the documentation for MITgcm's ACC channel simulation. At the date of writing, this documentation is also available online at https://mitgcm.readthedocs.io/en/latest/examples/reentrant_channel/reentrant_channel.html.

## `tools`
Miscellaneous helper functions.

## `verification`
Notebooks that inspect the MITgcm model output.

## `visualization`
For some of the figures in the manuscript. Others are located under the `ACC_particle_experiments` > `analysis` notebooks, or under `idealized_experiments/3D_3D_idealized_isopycnal_Markov1.ipynb`


