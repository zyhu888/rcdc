# **RCDC**: Robust Community Detection for noisy networks with Covariates

[![License:
GPL-3](https://img.shields.io/badge/License-GPL--3-green.svg)](LICENSE.md)

## RCDC

R package for Robust Community Detection for noisy networks with Covariates
(RCDC), integrating multiple noisy observations of a latent true network with
node-level covariates for community detection.

The main function is `RobustCDC()`, which fits RCDC over candidate community
numbers and selects the final model using WAIC.

## Installation

You can install `RCDC` from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("zyhu888/RCDC")
```

## Main Function

- `RobustCDC()` is the primary user-facing function. It fits candidate
  community numbers and returns the selected model.

## Parallel Runs

`RobustCDC()` evaluates candidate values of `K` serially. For large analyses,
run fixed-`K` jobs in parallel with `robustCDC_MCMC()` and combine the WAIC
results afterward. A SLURM/`parallel::mclapply()` template is included in the
parallel-runs vignette and at:

``` r
system.file("examples", "parallel-runs.R", package = "RCDC")
```

## References

References will be added after the manuscript is available.
