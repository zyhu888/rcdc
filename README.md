# **rcdc**: Robust Community Detection for noisy networks with Covariates

[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-green.svg)](LICENSE.md)

## RCDC

This R package performs robust community detection by jointly leveraging
multiple noisy observations of a latent true network and node-level covariates,
with the number of communities selected by minimizing WAIC.

## Installation

You can install `rcdc` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/rcdc")
```

## Main Functions

- `RobustCDC()` is the primary user-facing function. It fits candidate
  community numbers and returns the selected model.

- `robustCDC_MCMC()` runs the MCMC sampler for a fixed number of communities.

- `robustCDC_MCMC_saveallparameters()` runs the sampler and stores full
  parameter traces for detailed diagnostics.

- `dg()` simulates multiple noisy observations of a latent true network and
  covariates from user-specified parameters.

## Parallel Runs

`RobustCDC()` evaluates candidate values of `K` serially. For large analyses,
run fixed-`K` jobs in parallel with `robustCDC_MCMC()` and combine the WAIC
results afterward. A SLURM/`parallel::mclapply()` template is included in the
parallel-runs vignette and at:

``` r
system.file("examples", "parallel-runs.R", package = "rcdc")
```

## References

References will be added after the manuscript is available.
