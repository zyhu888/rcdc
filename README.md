RCDC: Robust Covariate-Assisted Community Detection
================

[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-green.svg)](LICENSE)

## Robust Covariate-Assisted Community Detection

The package implements Robust Covariate-Assisted Community Detection
(RCDC) for replicated noisy networks with node-level covariates.

This R package provides functions to:

- run the RCDC Markov chain Monte Carlo sampler;

- select the number of communities using WAIC;

- simulate network and covariate data for method evaluation;

- access supplemental helper routines used by the main fitting workflow.

## Installation

You can install `RCDC` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/RCDC")
```

## Main Functions

- `RobustCDC()` is the primary user-facing function. It fits candidate
  community numbers and returns the selected model.

- `robustCDC_MCMC()` runs the MCMC sampler for a fixed number of communities.

- `robustCDC_MCMC_saveallparameters()` runs the sampler and stores full
  parameter traces for detailed diagnostics.

- `sim_data_new_thetainput()` simulates replicated noisy networks and
  covariates from user-specified parameters.

## Parallel Runs

`RobustCDC()` evaluates candidate values of `K` serially. For large analyses,
run fixed-`K` jobs in parallel with `robustCDC_MCMC()` and combine the WAIC
results afterward. A SLURM/`parallel::mclapply()` template is included at:

``` r
system.file("examples", "parallel-runs.R", package = "RCDC")
```

## References

References will be added after the manuscript is available.
