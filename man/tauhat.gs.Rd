\name{taufuncs.gs}
\alias{tauhat.gs}
\title{ Internal Functions for Estimating tau }
\description{
These are internal functions used for calculating the scale parameter tau necessary for estimating the standard errors of coefficients for rank-regression.
}
\usage{
tauhat.gs(mids, counts, scores, bw = bw.nrd(mids))
}
\arguments{
  \item{mids}{vector. mid points of the interval}
  \item{counts}{vector. number of residuals in each bin}
  \item{scores}{vector. which scores were used}
  \item{bw}{bandwidth estimator.  default is rule of thumb by Scott (1992).}

}

\details{
Uses rule of thumb estimator in \code{bw.nrd} for speed.  Other estimators may be used if speed is less of a factor.  The default seems to work well based on simulations for large samples.
}

\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}

\author{
John Kloke, Joesph McKean, John Kapenga
Maintainer: John Kloke <kloke@biostat.wisc.edu>
}
