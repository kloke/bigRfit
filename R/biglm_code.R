#######################################################################
# Below modified from biglm 0.9-2.1 by Thomas Lumley June 2022 under GPL #
#######################################################################

"singcheck.bigqr" <-
function(bigQR){
  bigQR <- .Call("singcheckQR", bigQR)
}


"bigqr.init" <- function(p) {
  rval<-list(D=numeric(p), rbar=numeric(choose(p,2)),
             thetab=numeric(p),
             ss=0, checked=FALSE,
             tol=numeric(p))
  class(rval)<-"bigqr"
  rval
}

"coef.bigqr" <-
function(bigQR,nvar=NULL,...){
  p <- length(bigQR$D)
  if (is.null(nvar))
    nvar <- p
  if (nvar <1 | nvar >p) stop("Invalid value of `nvar'")

  if (!bigQR$checked)
    bigQR<-singcheck.bigqr(bigQR)
  
  tmp <- .Fortran("regcf", as.integer(p),
                  as.integer(p*p/2),
                  bigQR$D,
                  bigQR$rbar,
                  bigQR$thetab,
                  bigQR$tol,
                  beta=numeric(p),
                  nreq=as.integer(nvar),
                  ier=integer(1))

  if (tmp$ier!=0) stop("Error in REGCF: can't happen")

  tmp$beta
}


"update.bigqr" <-
function(bigQR, X, y, w=NULL,
                       singcheck=FALSE, add.intercept=FALSE){
  if (NCOL(X)+add.intercept!=length(bigQR$D))
    stop("Wrong number of columns")
  if (length(y)!=NROW(X))
    stop("Wrong number of rows")
  if (length(w)==0) w<-rep(1.0, length(y))
  if (length(y)!=length(w))
    stop("`weights' has wrong length")
  storage.mode(X)<-"double"
  storage.mode(y)<-"double"
  storage.mode(w)<-"double"
  bigQR<-.Call("updateQR",X, y, w, bigQR, add.intercept)
  
  if (singcheck)
    bigQR<-.Call("singcheckQR",bigQR);

  bigQR
}

