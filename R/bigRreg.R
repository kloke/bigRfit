bigRreg <- function(formula,data,...) { 

  call <- match.call()

  res <- bigRfit_xc(formula,data,intercept=TRUE,...) 

  res$call <- call

  res$betahat0 <- with(res, alphahat - drop(crossprod(xbar,betahat)))

  class(res) <- append(class(res),'bigRreg',0)

  res


}
