bigRfit <- function(x,y,B = min(max(length(y)^(2/3),1000),100000),scores=wscores,
  max.iter=100, eps=(.Machine$double.eps)^0.5) {

# note: normal scores don't work at the moment

if(length(y) <= 2000) stop("bigRfit requires at least 2000 records.  This is a job for rfit.")
# x is n x p design matrix
# y is n x 1 response vector
# B is the number of buckets

getScores.brf <- function(ehat,breaks,scores) {
  ngb <- hist(ehat,br=breaks,plot=FALSE)
  scores1 <- getScores(scores,c(0,cumsum(ngb$counts)/N))
  scoresvec <- (scores1[2:length(scores1)]+scores1[1:(length(scores1)-1)])/2
  cuts <- as.factor(cut(ehat,ngb$breaks,labels=FALSE))
  cutmat0 <- data.table(keys=levels(cuts),scores=scoresvec)
  cutsmat <- data.table(keys=cuts)
  scoremat <- merge(cutsmat,cutmat0,by='keys',all.x=TRUE,sort=FALSE)
  scoremat$scores
}

disp <- function( ehat, scores ) {
  sum( ehat*scores )
}

qrx <- qr(scale(x,center=TRUE,scale=FALSE))

ehat <- qr.resid(qrx,y)
breaks0 <- breaks <- quantile(ehat,seq(0,1,length=B))
breaks[1] <- -Inf
breaks[length(breaks)] <- Inf

scoresvec <- getScores.brf(ehat,breaks,scores=scores)
D0 <- disp(ehat,scoresvec)

i <- 0
while( i < max.iter ) {
  i <- i + 1
  
  tauhat1 <- tauhat.gs(breaks0,scores,bw.nrd(ehat))
  ehat <- ehat - tauhat1*qr.fitted(qrx,scoresvec)

  breaks0 <- breaks <- quantile(ehat,seq(0,1,length=B))
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  scoresvec <- getScores.brf(ehat,breaks,scores)

  D1 <- disp(ehat,scoresvec)

  if( (D0 - D1)/D0 < eps ) break()

  D0 <- D1

}
yhat <- y-ehat
betahat <- qr.coef(qrx,yhat)

alphahat <- median(ehat)

taustar <- taustar(ehat,qrx$rank)

list(coef=c(alphahat,betahat),fitted.values=yhat,residuals=ehat,tauhat=tauhat1,taustar=taustar,iter=i)

}
