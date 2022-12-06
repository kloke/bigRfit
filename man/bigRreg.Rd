\name{bigRreg}
\alias{bigRreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rank-based regression for big data.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bigRreg(formula, data, ...)
}
\arguments{
  \item{formula}{
formula
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

  \item{\dots}{
additional arguments passed to fitting routines
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Kloke, McKean (2023) "Nonparametric Statistical Methods using R"
}
\author{
John Kloke <johndkloke@gmail.com>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##
n <- 8000
x1 <- rt(n,9) ; x2 <- rt(n,9) 
y <- x1 + rt(n,9)
d <- data.frame(x1,x2,y)
bigRreg(y~1+x1+x2,data=d)  # fit model with intercept
bigRreg(y~x1+x2,data=d)  # same
bigRreg(y~x1+x2-1,data=d)  # same


## v0.8.1 December 2022 ##
## The function is currently defined as
function (formula, data, ...) 
{
    call <- match.call()
    res <- bigRfit_xc(formula, data, intercept = TRUE, ...)
    res$call <- call
    res$betahat0 <- with(res, alphahat - drop(crossprod(xbar, 
        betahat)))
    class(res) <- append(class(res), "bigRreg", 0)
    res
  }
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
