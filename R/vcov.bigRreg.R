vcov.bigRreg <- function(x) {

  pp1 <- length(x$xbar)+1

  K <- diag(pp1)
  K[1,(2:pp1)] <- -(x$xbar)

  K %*% vcov.bigRfit_xc1(x) %*% t(K)

}
