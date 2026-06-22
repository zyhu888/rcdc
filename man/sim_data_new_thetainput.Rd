\name{sim_data_new_thetainput}
\alias{sim_data_new_thetainput}
\title{Simulate RobustCDC Data}
\usage{
sim_data_new_thetainput(seed, n, K, M, p_cont, p_cat, F_s, theta_true,
  B, alpha, beta, mu, tau2, correct_rate)
}
\arguments{
\item{seed}{Random seed.}
\item{n}{Number of nodes.}
\item{K}{Number of communities.}
\item{M}{Number of replicated noisy networks.}
\item{p_cont}{Number of continuous covariates.}
\item{p_cat}{Number of categorical covariates.}
\item{F_s}{Number of categorical levels.}
\item{theta_true}{Node degree parameters.}
\item{B}{Community block matrix.}
\item{alpha}{False positive rate matrix.}
\item{beta}{False negative rate matrix.}
\item{mu}{Community means for continuous covariates.}
\item{tau2}{Covariance matrix for continuous covariates.}
\item{correct_rate}{Categorical covariate signal rates.}
}
\value{
A list containing noisy networks, covariates, and the data-generating truth.
}
\description{
Simulates replicated noisy networks and node covariates from user-specified parameters.
}
