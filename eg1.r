library(Rcpp)
library(data.table)

source('bigRfit.r')
source('tauhatB.r')
source('mytaustar.r')

N <- 10000
p <- 20
y <- rnorm(N)
x <- matrix(rnorm(N*p),ncol=p)

system.time(fit0<-rfit(y~x))
system.time(fit1<-bigRfit(x,y))

fit0$coef
fit1$coef

N <- 50000
p <- 40
y <- rnorm(N)
x <- matrix(rnorm(N*p),ncol=p)

system.time(fit0<-rfit(y~x))
system.time(fit1<-bigRfit(x,y))
