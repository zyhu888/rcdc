\name{robustCDC_MCMC_saveallparameters}
\alias{robustCDC_MCMC_saveallparameters}
\title{Run RobustCDC MCMC and Save Full Parameter Traces}
\usage{
robustCDC_MCMC_saveallparameters(A_tilde, X_cont = NULL, X_cat = NULL, K,
  n_iter = 3000, burn_in = 500, hyper, init = NULL)
}
\arguments{
\item{A_tilde}{A three-dimensional array of replicated noisy adjacency matrices.}
\item{X_cont}{Optional matrix of continuous node covariates.}
\item{X_cat}{Optional matrix of categorical node covariates.}
\item{K}{Number of communities.}
\item{n_iter}{Number of MCMC iterations.}
\item{burn_in}{Number of burn-in iterations.}
\item{hyper}{List of hyperparameters.}
\item{init}{Optional list of initial values.}
}
\value{
A list containing sampled networks, labels, parameters, burn-in, and log-likelihood matrices.
}
\description{
Runs the fixed-\code{K} sampler while retaining full parameter traces for diagnostics.
}
