% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_data.R
\name{sim.dat.fn}
\alias{sim.dat.fn}
\title{sim.dat.fn}
\usage{
sim.dat.fn(n, p, pi1, mu1, x.mu = 0, x.sd = 10, e.mu = 0, e.sd = 1)
}
\arguments{
\item{n}{Sample size.}

\item{p}{Number of covariates (including intercept).}

\item{pi1}{Proportion of outliers in the response variable.}

\item{mu1}{Expection of outliers.}

\item{x.mu}{Expection of design matrix X (default: 0).}

\item{x.sd}{Standard deviation of design matrix X (default: 10).}

\item{e.mu}{Expection of noise (default: 0).}

\item{e.sd}{Standard deviation of noise (default: 1).}
}
\value{
Simulated data and its parameters.
}
\description{
To simulation linear regression dataset
}
\author{
Meng Wang
\email{mengw1@stanford.edu}
}
