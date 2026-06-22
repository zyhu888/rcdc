\name{robustCDC_MCMC}
\alias{robustCDC_MCMC}
\title{Run the RobustCDC MCMC Sampler}
\usage{
robustCDC_MCMC(A_tilde, X_cont = NULL, X_cat = NULL, K,
  n_iter = 3000, burn_in = 500, hyper, init = NULL)
}
\arguments{
\item{A_tilde}{A three-dimensional array of noisy adjacency matrices observed from the same latent true network.}
\item{X_cont}{Optional matrix of continuous node covariates.}
\item{X_cat}{Optional matrix of categorical node covariates.}
\item{K}{Number of communities.}
\item{n_iter}{Number of MCMC iterations.}
\item{burn_in}{Number of burn-in iterations.}
\item{hyper}{List of hyperparameters.}
\item{init}{Optional list of initial values.}
}
\value{
A list containing sampled labels, burn-in, and pointwise log-likelihood matrices.
}
\description{
Runs the fixed-\code{K} RobustCDC Markov chain Monte Carlo sampler.
}
