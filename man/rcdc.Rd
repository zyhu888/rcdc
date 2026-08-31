\name{rcdc}
\alias{rcdc}
\title{Fit RCDC and Select the Number of Communities}
\usage{
rcdc(A_tilde, X_cont = NULL, X_cat = NULL, K = 2:10,
  n_iter = 3000, burn_in = 500, hyper, init = NULL, verbose = TRUE)
}
\arguments{
\item{A_tilde}{A three-dimensional array of noisy adjacency matrices observed from the same latent true network.}
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
An object of class \code{rcdc} with selected \code{K}, WAIC, labels, and a WAIC table.
}
\description{
Runs \code{rcdc_MCMC()} over candidate values of \code{K}, computes WAIC, and returns the model with the smallest WAIC.
}
\details{
\code{rcdc()} is the primary interface for fitting RCDC. For each candidate value in \code{K},
the function calls \code{rcdc_MCMC()}, stores the pointwise log-likelihood
after burn-in, computes WAIC, and keeps the clustering estimate from the model
with the smallest WAIC.

The input \code{A_tilde} should be an \code{n} by \code{n} by \code{M} array,
where each slice is a noisy adjacency matrix on the same set of
nodes. Continuous covariates are supplied through \code{X_cont}; categorical
covariates are supplied through \code{X_cat}. The hyperparameter list controls
the priors for the latent network, observation error rates, covariate
distributions, and community proportions.

The returned label vector \code{z} is obtained from the post-burn-in MCMC
samples by taking the posterior modal community label for each node.
}
\examples{
library(rcdc)
library(aricode)

seed <- 123
set.seed(seed)

n <- 60
K0 <- 3
M <- 4
burn <- 200
niter <- 500

rho <- c(0.10, 0.30)
zeta <- c(1.0, 1.0)

hyper <- list(
  sigma = 10,
  a_alpha = 1,
  b_alpha = 1,
  a_beta = 1,
  b_beta = 1,
  mu0 = rep(0, K0),
  iota0 = rep(1, K0),
  a_tau = rep(1, K0),
  b_tau = rep(0.01, K0),
  F_cat = c(K0, K0),
  gamma0 = c(1, 1),
  d = 1
)

B <- matrix(-1.5, K0, K0)
diag(B) <- 0

mu <- diag(K0) - matrix(1 / K0, K0, K0)
mu <- mu * (10 / sqrt(2))

tau2 <- diag(3^2, K0)
theta <- rnorm(n, mean = zeta[1], sd = zeta[2])
alpha <- matrix(rho[1], K0, K0)
beta <- matrix(rho[2], K0, K0)
cat_rate <- c(0.7, 0.7)

sim <- dg(
  seed = seed,
  n = n,
  K = K0,
  M = M,
  p_cont = K0,
  p_cat = 2,
  F_s = K0,
  theta_true = theta,
  B = B,
  alpha = alpha,
  beta = beta,
  mu = mu,
  tau2 = tau2,
  correct_rate = cat_rate
)

fit <- rcdc(
  A_tilde = sim$A_tilde,
  X_cont = sim$X_cont,
  X_cat = sim$X_cat,
  K = 2:5,
  n_iter = niter,
  burn_in = burn,
  hyper = hyper
)

print(fit)

z_hat <- fit$z
z_true <- sim$truth$z

ARI(z_hat, z_true)
NMI(z_hat, z_true)
length(unique(z_hat))
fit$waic_table
}
\references{
XXX. References will be added after the manuscript is available.
}
