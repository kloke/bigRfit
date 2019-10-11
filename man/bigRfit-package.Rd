\name{bigRfit-package}
\alias{bigRfit-package}
\docType{package}
\title{bigRfit - rank-based estimation for linear models with big data}
\description{
Rank-based (R) estimation for linear models with big data.
Requires at least 200 records (for smaller datasets use Rfit).
}
\details{
\tabular{ll}{
Package: \tab bigRfit\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2016-06-28\cr
License: \tab GPL (>=2)\cr
}
}
\author{
John Kloke, Joesph McKean, John Kapenga
Maintainer: John Kloke <johndkloke@gmail.com>
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\keyword{ package }
\seealso{
code{\link[Rfit]{Rfit}} ~~
}
\examples{
n <- 10^6; p <- 10
x <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
#system.time(fit1 <- bigRfit(x,y))
#summary(fit1)
}
