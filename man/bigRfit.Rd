\name{bigRfit}
\alias{bigRfit}
\title{Rank-based estimation for linear models with big data.}
\description{
An implementation of rank-based estimation when working with big data.
Uses step scores to increase the speed.
Requires at least 2003 records in the dataset.
}
\usage{
bigRfit(x, y, B = 1001, scores = Rfit::wscores, max.iter = 100, eps = (.Machine$double.eps)^0.625)
}
\arguments{
  \item{x}{n by p design matrix }
  \item{y}{n by 1 response vector }
  \item{B}{number of breaks (number of 'buckets' + 1)}
  \item{scores}{an object of class scores}
  \item{max.iter}{maximum number of iterations}
  \item{eps}{stopping criteria}
}
\value{
\item{coefficients}{estimated regression coefficents with intercept}
\item{residuals}{the residuals, i.e. y-yhat}
:\item{fitted.values}{ yhat = x betahat}
\item{scores}{ score function used in estimation}
\item{x}{design matrix w/ intercept added (ie cbind(1,x))}
\item{y}{original response vector}
\item{tauhat}{estimated value of the scale parameter tau}
\item{taushat}{estimated value of the scale parameter tau_s}
\item{symmetric}{currently set to TRUE (needed for helper functions down the line)}
\item{iter}{number of iterations}
\item{D1}{final model dispersion}
\item{D0}{null model dispersion}
\item{converage}{convergance status (logical)}
\item{qrx1}{result of call to qr of cbind(1,x)}
}
\references{
Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
}
\author{John Kloke \email{johndkloke@gmail.com}}

\seealso{\code{\link{Rfit}} }
\examples{
n <- 10^4; p <- 10
x <- matrix(rnorm(n*p),ncol=p)
y <- rnorm(n)
#system.time(fit0 <- rfit(y~x))
#summary(fit0)
#system.time(fit1 <- bigRfit(x,y))
#summary(fit1)
}
\keyword{nonparametric}
\keyword{robust}
\keyword{regression}

