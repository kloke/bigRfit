\name{bigRfit_biglm_xc}
\alias{bigRfit_biglm_xc}
\alias{bigRfit_xc}
\title{
Big data rank-based estimation based on centered design matrix.
}
\description{
Fits the regression coeficients for the linear regression model with a centered design matrix.
Algorithm designed for big data sets which uses a partial ranking and step scores.
}
\usage{
bigRfit_biglm_xc(formula, data, intercept = FALSE, yhat0 = NULL, ehat0 = NULL, B = 1000, scores = wscores, max.iter = 50, eps = sqrt(.Machine$double.eps), TAU = "DT", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class formula
}
  \item{data}{
data frame
}
  \item{intercept}{
indicator to request estimate of alpha (FALSE by default)
}
  \item{yhat0}{
optional initial estimate of responses
}
  \item{ehat0}{
optional initial estimate of residuals
}
  \item{B}{
number of bins to use (default of 1000)
}
  \item{scores}{
an object of class 'scores' (wscores by default)
}
  \item{max.iter}{
maximum number of iterations
}
  \item{eps}{
specify tol
}
  \item{TAU}{
version of estimation routine for scale parameter.  DT for approximate binned estimate (using data.table), F0 for Fortran using full set of residuals, N for none
}
  \item{\dots}{
Additional arguments.
}
}
\details{
Given a formula y ~ x fits the model y ~ xc 
where xc is the centered x.  
Estimation of an intercept is not needed.  
However, an estimate of alpha in the model
y = alpha 1 + xc beta + e
may be requested using the option intercept=TRUE.

Rank regression estimates beta by minimizing the dispersion functions:
D(beta) = a(R(y-xc beta))^T (y-xc beta)
where R denotes rank and a is a non-decreasing score function.

In this implementation,
instead of using a full ranking, the observations are binned into 
bins of approximate equal size. 
The parameter B specifies the number of bins to use.
B must be <= the number of observations.
Step scores are used so that observations within the same bin are assigned the same score.

This function is (mostly) intended to be an internal function.
See bigRreg for general purpose regression. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
John Kloke <johndkloke@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bigRreg}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
