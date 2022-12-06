bigRfit_xc <- function(...) bigRfit_biglm_xc(...)

bigRfit_biglm_xc <- function(formula, data, intercept=FALSE, yhat0=NULL, ehat0=NULL, B=1000, scores=wscores,
                          max.iter = 50, eps = sqrt(.Machine$double.eps), TAU='DT',...
                         ) {

# returns estimated coef for the centered design model
# i.e., fits intercept and regression coef in 
#   y = alpha * one + Xc * beta + error
# where Xc is the centered version of the provided design matrix 
# centering by the function scale with arguments scale=FALSE, center=TRUE

# bandwidth approximated using midpoints of bins

if( (!intercept) & (TAU == 'F0') ) stop("Estimate of alpha required for F0 estimate of TAU.")

call<-match.call()

# set max.iter to NULL to have it choose
if(is.null(max.iter)) max.iter <- 25*(n<10^5) + 20*(n<10^7) + 5

############################
### function definitions ###
############################
bw_bigRfit0 <- function(x,n) {
# modified version of bw.nrd under GPL

  r <- quantile(x, c(0.25, 0.75))
  h <- (r[2L] - r[1L])/1.34
  1.06 * min(sqrt(var(x)), h) * n^(-1/5)

}

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


get_breaks <- function(e,B1,EPS=0.001) {
  breaks <- unique(quantile(e,seq(0,1,length=B1)))
  nbs <- length(breaks)
  ind <- c(1,nbs)
  s <- EPS*max( sqrt(.Machine$double.eps), diff(breaks[c(0.25,0.75)*nbs]))
  breaks[ind] <- breaks[ind] + EPS*diff(breaks[c(0.25,0.75)*nbs])*c(-1,1)  # one hundredth of one percent
  breaks <- unique(breaks)
  breaks
}

get_mids <- function(breaks) { 
  nb <- length(breaks)
  0.5*(breaks[1:(nb-1)]+breaks[2:nb])
}

get_scoreDT <- function(ehatDT) {

  scoreMat <- ehatDT[,.(count=.N),by=bin]
  setkey(scoreMat,bin)  #sorts by bin

  scoreMat[,rnkH := cumsum(count)]
  scoreMat[,rnkL := rnkH-count+1]
  scoreMat[,rnkH := rnkH/(n+1)]
  scoreMat[,rnkL := rnkL/(n+1)]

  scoreMat[,scrs := 0.5*(getScoresNS(scores,rnkH)+getScoresNS(scores,rnkL))]
  m <- scoreMat[,sum(count*scrs)/n]
  scoreMat[,scrs := scrs - m]
  sd <- sqrt(scoreMat[,sum(scrs*scrs*count/(n+1)) ])
  scoreMat[,scrs := scrs / sd]
  scoreMat[,scrsD := 0.5*(getScoresDerivNS(scores,rnkH)+getScoresDerivNS(scores,rnkL))/sd]

  return(data.table(scoreMat))

}

dispDT <- function(scoreMat,ehatDT) {
  merge(scoreMat,ehatDT[,.(D=sum(ehat)),bin],by='bin',sort=FALSE)[,sum(D*scrs)]
}

tauhatDT <- function(scoreMat,mids,n,h) { 
#  h <- bw_bigRfit(mids,n,scoreMat[,count])
  pc <- scoreMat[,count]/n
  D <- (dnorm(outer(mids,mids,'-')/h)/h)%*%pc
  tauhat0 <- 1/drop(crossprod(pc*unlist(scoreMat[,'scrsD']),D))
  return(tauhat0)
}


############################
############################


#######################################################################
# Below modified from biglm 0.9-2.1 by Thomas Lumley June 2022 under GPL #
off <- 0
w<-NULL
ttxc<-terms(update(formula, ~ . -1))
mfxc<-model.frame(ttxc,data)
#x<-model.matrix(tt,mf)
#x<-as.matrix(x[, colnames(x) != "(Intercept)"])
#x1 <- as.matrix(cbind(rep(1, nrow(x)), x))
#xc<-scale(x,scale=FALSE,center=TRUE)
xc<-scale(model.matrix(ttxc,mfxc),scale=FALSE,center=TRUE)

y <- model.response(mfxc)
scale.y <- FALSE
if(scale.y) {
  scl <- 10^floor(log10(diff(range(y))))
  y <- y/scl
}

qrx0<-bigqr.init(NCOL(xc))
if(is.null(ehat0)) { 
  if(!is.null(yhat0)) {
    ehat0 <- y-yhat0
  } else {
    qr0<-update(qrx0,xc,y-off,w)
    betahat0 <- coef.bigqr(qr0)
    ehat0 <- y-drop(xc%*%betahat0)
  }
}
#######################################################################
n <- nrow(xc)
p <- ncol(xc)

if(is.null(max.iter)) max.iter <- 20*(n<10^4) + 15*(n<10^5) + 10*(n<10^6) + (n<10^7) + 1

update2disp <- function(ehat,B1) {
  breaks <- get_breaks(ehat,B1)
  mids <- get_mids(breaks)
  ehatDT <- data.table(ehat=ehat)
  ehatDT[,bin:=cut(ehat,breaks,label=FALSE)]
  scoreMat <- get_scoreDT(ehatDT)
  D <- dispDT(scoreMat,ehatDT)
  return(list(mids=mids,ehatDT=data.table(ehatDT),scoreMat=scoreMat,D=D))
}

#########################################
### values based on initial residuals ###
#########################################
breaks <- get_breaks(ehat0,B+1)
mids <- get_mids(breaks)

# set up data table for resids
ehatDT <- data.table(ehat=ehat0)
ehatDT[,bin:=cut(ehat,breaks,label=FALSE)] 

scoreMat <- get_scoreDT(ehatDT)
D0 <- dispDT(scoreMat,ehatDT)
tauhat <- tauhatDT(scoreMat,mids,n,h=bw_bigRfit(mids,n,scoreMat[,count]))
#########################################
#########################################

converge <- FALSE

for( iter in seq_len(max.iter) ) {  ## start main loop ##

############################
### Newton-step estimate ###
############################
#aVec <- unlist(merge(ehatDT,scoreMat[,c('bin','scrs')],by='bin',sort=FALSE)[,'scrs'],use.name=FALSE)
#qr1<-update(qrx0,xc,aVec,w)
qr1<-update(qrx0,xc,
            unlist(merge(ehatDT,scoreMat[,c('bin','scrs')],by='bin',sort=FALSE)[,'scrs'],use.name=FALSE),
            w)
rm(ehatDT)
dir1 <- tauhat*coef(qr1)

betahat1 <- betahat0 + dir1

foo1 <- update2disp(y-drop(xc%*%betahat1),B+1L)
D1 <- foo1$D
#########################
#########################

###################
### backtracking ##
###################
backtrack.denom <- crossprod(dir1)/tauhat

backtrack.c <- 0.0001
backtrack.alpha <- 1
backtrack.rho <- 0.5

repeat{

  betahat2 <- betahat1 + backtrack.alpha*dir1
  foo2 <- update2disp(y-drop(xc%*%betahat2),B+1L)
  D2 <- foo2$D

  if( (D2-D1)/backtrack.denom <= backtrack.c*backtrack.alpha) (break)()
  backtrack.alpha <- backtrack.alpha*backtrack.rho

}

if( D1 < D2 ) {
cat("Note: setting estimate to previous backtrack.\n")
  D2 <- D1
  betahat2 <- betahat1
  foo2 <- foo1
}

###################
###################

###########################
## check for convergence ##
###########################

if( (D0 - D2)/D0 < eps ) {
  converge <- TRUE
  break()
}

D0 <- D2
betahat0 <- betahat2
ehatDT <- foo2$ehatDT

} ## end main loop ##

D1 <- D2
betahat <- betahat2
scoreMat <-foo2$scoreMat
mids <- foo2$mids

# estimate of tau
approximate.final.h <- TRUE
if(approximate.final.h) {
  final.h <- bw_bigRfit(mids,n,scoreMat[,count])
} else {
  final.h <- bw_bigRfit0(ehatDT[,ehat],n)
}

xnames <- colnames(xc)

coef <- NULL
alphahat <- taushat <- NULL
if( intercept ) {
xnames <- append(xnames,'Intercept',0)
y1 <- drop(xc%*%betahat)
ehat <- y - y1
alphahat <- median(ehat)
ehat <- ehat-alphahat
taushat <- taustar(ehat,p+1)
}


if(TAU == 'N' ) {
  tauhat <- NULL
} else { 
  if(TAU == 'DT') {
    tauhat <- tauhatDT(scoreMat,mids,n,h=final.h)
  } else {
    tauhat <- gettauBDF0(scoreMat,ehat,n=n,p=p)
  }
}


## null model disp ##
D0 <- update2disp(y,B+1)$D


res <- list(coef=coef,alphahat=alphahat, betahat=betahat, D1=D1,D0=D0,tauhat=tauhat, taushat=taushat, qrxc=qr1,n=n,iter=iter,converge=converge, B=B, n=n, xbar=attr(xc,'scaled:center'),TAU=TAU,xnames=xnames)
res$call <- call
class(res) <- list('bigRfit_xc','bigRfit')

return(res)

}
