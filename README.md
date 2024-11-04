# Analysis of sea otter foraging observations in Glacier Bay, Alaska

This repository contains code to investigate spatio-temporal patterns in sea otter prey selection over the course of sea otter expansion across Glacier Bay from 1993 to 2019. 

All code is stored in the `src` directory and includes:
-  `bout_level_nb.stan`: file that defines the model in Stan
-  `run_stan.r`: script that loads data and runs the MCMC using `rstan`
-  `process_data.r`: script to process the foraging data and align with estimates of local sea otter abundance
-  `plots.r`: script to generate all the plots in the manuscript from output of `run_stan.r`
-  `eda.r`: script to generate exploratory figures and summaries

The scripts above write and read MCMC output to an `output` directory.

Data and model inputs are expected to be stored in a `data` directory.
-  Sea otter foraging observations: https://doi.org/10.5066/P9T6UB5G
-  Model-based estimates of sea otter abundance: https://doi.org/10.5061/dryad.vt4b8gtx6
