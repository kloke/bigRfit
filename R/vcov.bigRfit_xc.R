vcov.bigRfit_xc0 <- function(x) {

  p <- length(x$qr$D)
  R <- diag(p)
  R[row(R) > col(R)] <- x$qr$rbar
  R <- t(R)
  R <- sqrt(x$qr$D)*R
  ok <- x$qr$D!=0
  R[ok,ok]<-chol2inv(R[ok,ok])
  R[!ok,]<-R[,!ok]<-NA

  R*(x$tauhat*x$tauhat)

}

vcov.bigRfit_xc1 <- function(x) {

if(is.null(x$alphahat)) warning('No estimate of alpha found.')

  pp1 <- length(x$xbar) + 1

  V <- matrix(0,nrow=pp1,ncol=pp1)
  V[1,1] <- (x$taushat^2)/x$n
  V[2:pp1,2:pp1] <- vcov.bigRfit_xc0(x)

  V

}

vcov.bigRfit_xc <- function(x) {
  vcov.bigRfit_xc0(x)
#  if(is.null(x$alphahat)) return(vcov.bigRfit_xc0(x))
#  else return(vcov.bigRfit_xc1(x))
}
