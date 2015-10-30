bigRfit <- function(x,y,B = min(max(length(y)^(2/3),1000),100000),
  max.iter=100, eps=(.Machine$double.eps)^0.5) {

# x is n x p design matrix
# y is n x 1 response vector
# B is the number of buckets

wilscrs <- function(x) sqrt(12)*(x-0.5)

getScores <- function(ehat,breaks) {
  ngb <- hist(ehat,br=breaks,plot=FALSE)
  scores1 <- wilscrs(c(0,cumsum(ngb$counts)/N))
  scores <- (scores1[2:length(scores1)]+scores1[1:(length(scores1)-1)])/2
  cuts <- as.factor(cut(ehat,ngb$breaks,labels=FALSE))
  cutmat0 <- data.table(keys=levels(cuts),scores)
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

scores <- getScores(ehat,breaks)
D0 <- disp(ehat,scores)

i <- 0
while( i < max.iter ) {
  i <- i + 1
  
  tauhat1 <- tauhat(breaks0,bw.nrd(ehat))
  ehat <- ehat - tauhat1*qr.fitted(qrx,scores)

  breaks0 <- breaks <- quantile(ehat,seq(0,1,length=B))
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  scores <- getScores(ehat,breaks)

  D1 <- disp(ehat,scores)

  if( (D0 - D1)/D0 < eps ) break()

  D0 <- D1

}
yhat <- y-ehat
betahat <- qr.coef(qrx,yhat)

alphahat <- median(ehat)

taustar <- taustar(ehat,qrx$rank)

list(coef=c(alphahat,betahat),fitted.values=yhat,residuals=ehat,tauhat=tauhat,taustar=taustar,iter=i)

}
