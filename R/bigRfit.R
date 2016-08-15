bigRfit <-
function(x,y,B = 101, scores=Rfit::wscores,
  max.iter=100, eps=(.Machine$double.eps)^0.5) {

# should scale y by sd 

# note: normal scores don't work at the moment

if(length(y) <= 200) stop("bigRfit requires at least 200 records.  This is a job for rfit.")
# x is n x p design matrix
# y is n x 1 response vector
# B is the number of buckets (plus 1?)

getScores.brf <- function(ehat,breaks,scores) {
#  breaks <- unique(breaks)
  ngb <- hist(ehat,br=breaks,plot=FALSE)
  scores1 <- getScores(scores,c(0,cumsum(ngb$counts)/sum(ngb$counts)))
  scoresvec <- (scores1[2:length(scores1)]+scores1[1:(length(scores1)-1)])/2
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
  sum( ehat*scores )
}

get_breaks <- function(ehat,B) {
  breaks <- quantile(ehat,seq(0,1,length=B))
  ind <- c(1,length(breaks))
  breaks[ind] <- breaks[ind] + eps*sd(breaks)*c(-1,1)
  unique(breaks)
}

qrx <- qr(scale(x,center=TRUE,scale=FALSE))

ehat <- qr.resid(qrx,y)

breaks <- get_breaks(ehat,B)
scrs <- getScores.brf(ehat,breaks,scores=scores)

#scoresvec <- getScores.brf(ehat,breaks,scores=scores)
D0 <- disp(ehat,scrs$scorevec)

converge <- FALSE
i <- 0
while( i < max.iter ) {
  i <- i + 1
  
  tauhat1 <- tauhat.gs(scrs$mids,scrs$counts,scores,bw.nrd(ehat))
  dir1 <- -tauhat1*qr.fitted(qrx,scrs$scorevec)
  ehat1 <- ehat + dir1

  breaks <- get_breaks(ehat1,B)
#  breaks[c(1,length(breaks))] <- c(-Inf,Inf)
#  scoresvec <- getScores.brf(ehat,breaks,scores)
  scrs <- getScores.brf(ehat1,breaks,scores=scores)

  D1 <- disp(ehat1,scrs$scorevec)

  if( D1 > D0 ) {
    ehat1 <- ehat + 0.5 * dir1
  } else {

    if( (D0 - D1)/D0 < eps ) { 
      converge <- TRUE
      break()
    }

    D0 <- D1

  }

  ehat <- ehat1


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

res <- list(coefficients=fit.ls$coefficients,
  fitted.values=yhat,residuals=ehat,x=cbind(1,x),y=y,
  tauhat=tauhat1,taushat=taustar,
  qrx1=fit.ls$qr, D1=D1, D0=D0, scores=scores, symmetric=TRUE,
  iter=i,converge=converge)

class(res) <- list("bigRfit","rfit")

res

}
