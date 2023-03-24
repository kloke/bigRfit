bw_bigRfit <- function(x,n,nvec) {

# modified version of bw.nrd under GPL
# using binned data to approximate the bandwidth
# with a large number of bins and large samples
#   intuitively it should provide a good approximation
# x - midpoints of bins
# n - number of observations
# nvec - counts per bin

# approximation based on already approximate quantiles
  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2L] - r[1L])/1.34
#  1.06 * min(sqrt(var(x)), h) * n^(-1/5)

# approximation based on the binned data
  xbar <- drop(crossprod(x,nvec))/n
  v <- sum(nvec*(x-xbar)^2) / (n-1)

  1.06 * min(sqrt(v), h) * n^(-1/5)
}

approxtauDT <- function(scoreMat,mids,n,h=bw_bigRfit(mids,n,scoreMat[,count])) { 
#  h <- bw_bigRfit(mids,n,scoreMat[,count])
  pc <- scoreMat[,count]/n
  D <- (dnorm(outer(mids,mids,'-')/h)/h)%*%pc
  tauhat0 <- 1/drop(crossprod(pc*unlist(scoreMat[,'scrsD']),D))
  return(tauhat0)
}

approxtauAH <- approxtauDT
