bigRfit <-
function(x,y,B = 1001, scores=Rfit::wscores,
  max.iter=100, eps=(.Machine$double.eps)^0.625) {

# should scale y by sd 

# input checks
if(length(y) <= 2002) stop("bigRfit requires at least 2002 records.  This is a job for rfit.")
if( length(y)/B < 2 ) stop("Number of bins requested too large. Consider using rfit or reducing number of bins.")

# potiential arguments to function
trace <- 1
yhat0 <- y
#qrx <- qr(scale(x,center=TRUE,scale=FALSE))

if(length(y) != nrow(x)) stop("response vector (y) not in the same dimension as design matrix (x)")
#n <- length(y)

# x is n x p design matrix
# y is n x 1 response vector
# B is the number of bins + 1

getScores.brf <- function(ehat,breaks,scores) {
#  breaks <- unique(breaks)
  ngb <- hist(ehat,br=breaks,plot=FALSE)
  rank1 <- cumsum(ngb$counts)      #high rank
  rank0 <- rank1 - ngb$counts + 1  #low rank
  ave_rank <- (rank1 + rank0)/2
  scoresvec <- getScores(scores,ave_rank/(length(ehat)+1))
#  scoresvec <- (scores1[2:length(scores1)]+scores1[1:(length(scores1)-1)])/2
#  cuts <- as.factor(cut(ehat,ngb$breaks,labels=FALSE))
  cuts <- cut(ehat,ngb$breaks,labels=FALSE)
  lc <- sort(unique(cuts))
  scoresvec <- scoresvec[lc]
#  cutmat0 <- data.table(keys=lc,scores=scoresvec)
#  cutsmat <- data.table(keys=as.numeric(cuts))
  cutmat0 <- data.table(keys=lc,scores=scoresvec)
  cutsmat <- data.table(keys=cuts)
  scoremat <- merge(cutsmat,cutmat0,by='keys',all.x=TRUE,sort=FALSE)
  scorevec <- scoremat$scores - mean(scoremat$scores)
  list(scorevec=scorevec,mids=ngb$mids,counts=ngb$counts)
}

disp <- function( ehat, scores ) {
  drop(crossprod(ehat,scores))
#  sum( ehat*scores )
}

get_breaks <- function(ehat,B) {
  breaks <- quantile(ehat,seq(0,1,length=B))
  ind <- c(1,length(breaks))
  breaks[ind] <- breaks[ind] + eps*sd(breaks)*c(-1,1)
  unique(breaks)
}

qrx <- qr(scale(x,center=TRUE,scale=FALSE))

ehat <- qr.resid(qrx,yhat0)
# rm(yhat0)

breaks <- get_breaks(ehat,B)
scrs <- getScores.brf(ehat,breaks,scores=scores)

#scoresvec <- getScores.brf(ehat,breaks,scores=scores)
D1 <- disp(ehat,scrs$scorevec)

converge <- FALSE
i <- 0

ehat0 <- ehat
D0 <- D1

while( i < max.iter ) {
  i <- i + 1
 if( trace ) cat(paste0("iter = ", i, "\n"))
  tauhat1 <- tauhat.gs(scrs$mids,scrs$counts,scores,bw.nrd(ehat))
if( trace ) cat(paste0("  tauhat = ", tauhat1, "\n"))
  dir1 <- -tauhat1*qr.fitted(qrx,scrs$scorevec)

#############
# backtrack #
#############
backtrack.denom <- crossprod(dir1)/tauhat1

backtrack.c <- 0.0001
backtrack.alpha <- 1
backtrack.rho <- 0.5
repeat {

if(trace) cat(paste0('Backtrack : alpha = ', backtrack.alpha, "\n"))
  ehat <- ehat0 + backtrack.alpha*dir1
  breaks <- get_breaks(ehat,B)
#  breaks[c(1,length(breaks))] <- c(-Inf,Inf)
#  scoresvec <- getScores.brf(ehat,breaks,scores)
  scrs <- getScores.brf(ehat,breaks,scores=scores)

  D1 <- disp(ehat,scrs$scorevec)

if(trace) cat(paste0('  D1 = ', D1, "\n"))
if(trace) cat(paste0('  D0 = ', D0, "\n"))
  if( (D1 - D0)/backtrack.denom <= backtrack.c*backtrack.alpha) (break)()
  backtrack.alpha <- backtrack.alpha*backtrack.rho

}
#############
#############
    if( (D0 - D1)/D0 < eps ) { 
      converge <- TRUE
      break()
    }

  D0 <- D1
  ehat0 <- ehat
}

if(!converge) cat('Note: convergance criteria not met.  iter = ', i, "\n")

if( D0 < D1 ) {
cat("Note: setting estimate to previous iteration.\n")
  D1 <- D0
  ehat <- ehat0
}

alphahat <- median(ehat)
ehat <- ehat - alphahat
yhat <- y-ehat
fit.ls <- lsfit(x,yhat,intercept=TRUE)

taustar <- taustar(ehat,qrx$rank)

  breaks <- quantile(y,seq(0,1,length=B))
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  breaks <- unique(breaks)
  scrs <- getScores.brf(y,breaks,scores=scores)
#  scoresvec <- getScores.brf(y,breaks,scores)
  D0 <- disp(y,scrs$scorevec)


# Huber's DF correction
dftauhat <- function(tauhat,ehat,param,p){

        epshc <- 0.000001
        n <- length(ehat)
#        ind <- rep(0,n)
#        ind[abs(ehat/tauhat) < param] <- 1
#        hubcor <- sum(ind)/n
#
  hubcor <- mean( abs(ehat/tauhat) < param )
        if(hubcor < epshc){hubcor <- epshc}
        corr <- sqrt(n/(n-p))*(1 + ((p/n)*((1-hubcor)/hubcor)))
        tauhat <- tauhat*corr
        tauhat
}

huber <- FALSE
if(huber) tauhat1 <- dftauhat(tauhat1,ehat,param=2,qrx$rank)

res <- list(coefficients=fit.ls$coefficients,
  fitted.values=yhat,residuals=ehat,x=cbind(1,x),y=y,
  tauhat=tauhat1,taushat=taustar,
  qrx1=fit.ls$qr, D1=D1, D0=D0, scores=scores, symmetric=TRUE,
  iter=i,converge=converge)

class(res) <- list("bigRfit","rfit")

res

}
