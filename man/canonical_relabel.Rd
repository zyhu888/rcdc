\name{canonical_relabel}
\alias{canonical_relabel}
\title{Canonical Relabeling of Community Labels}
\usage{
canonical_relabel(z, K)
}
\arguments{
\item{z}{Integer vector of community labels.}
\item{K}{Number of communities.}
}
\value{
A list with relabeled labels \code{z} and the permutation order \code{ord}.
}
\description{
Relabels communities by first appearance to reduce label switching.
}
