\name{RobustCDC}
\alias{RobustCDC}
\title{Fit RobustCDC and Select the Number of Communities}
\usage{
RobustCDC(A_tilde, X_cont = NULL, X_cat = NULL, K = 2:10,
  n_iter = 3000, burn_in = 500, hyper, init = NULL, verbose = TRUE)
}
\arguments{
\item{A_tilde}{A three-dimensional array of replicated noisy adjacency matrices.}
\item{X_cont}{Optional matrix of continuous node covariates.}
\item{X_cat}{Optional matrix of categorical node covariates.}
\item{K}{Candidate numbers of communities.}
\item{n_iter}{Number of MCMC iterations.}
\item{burn_in}{Number of burn-in iterations.}
\item{hyper}{List of hyperparameters.}
\item{init}{Optional list of initial values.}
\item{verbose}{Logical; print progress messages.}
}
\value{
An object of class \code{robustCDC} with selected \code{K}, WAIC, labels, and a WAIC table.
}
\description{
Runs \code{robustCDC_MCMC()} over candidate values of \code{K}, computes WAIC, and returns the model with the smallest WAIC.
}
