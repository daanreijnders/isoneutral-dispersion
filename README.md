# Simulating Lagrangian Subgrid-Scale Dispersion on Neutral Surfaces in the Ocean
This repository holds the data generation and analysis scripts for the manuscript 'Simulating Lagrangian Subgrid-Scale Dispersion on Neutral Surfaces in the Ocean'. 

## `ACC_mitgcm_config`
Contains MITgcm configuration files for the coarse and fine model, including spinups. This includes the compilation code and namelists.

## `ACC_particle_experiments`
Contains Python scripts for the Lagrangian simulations carried out with [Parcels](https://github.com/OceanParcels/parcels) (experiments with a tracer patch, and with particles initialized on a lattice), and analysis scripts pertaining to these simulations and their output.

## `data_processing`
Contains scripts and notebooks for processing the data. This includes creating coarsened fields, computing discrete derivatives, and computing the quantities necessary for the LS parameterization.

## `idealized_experiments`
Experiments using the idealized set-up. 

## `kernels`
Functions used by [Parcels](https://github.com/OceanParcels/parcels) to drive the Markov models.

## `tools`
Miscellaneous helper functions.

## `verification`
Notebooks that inspect the MITgcm model output.

## `visualization`
For some of the figures in the manuscript. Others are located under the `ACC_particle_experiments` > `analysis` notebooks.

