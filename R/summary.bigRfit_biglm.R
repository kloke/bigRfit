summary.bigRfit <- function(object,...) {

# object of class bigRfit_xc or bigRreg 
#   i.e., expects object$coef to be NULL (exits with error)

  if(!is.null(object$coef)) stop("Summary not yet implemented for bigRfit.")

  pp1 <- length(object$betahat) + 1

  se <-  sqrt(diag(vcov(object)))

  if(is.null(object$betahat0)) {
    est <- object$betahat
  } else {
    est <- with(object,c(betahat0,betahat))
  }

  if( length(est) != length(se) ) stop("Estimate and Standard Error not equal dim.")

#  se <-  sqrt(diag(vcov(object)))
  tstat <- est/se

  pval <- 2 * pt(-abs(tstat), n - pp1)
#    coef <- cbind(est, ses, tstat, pval)
#    colnames(coef) <- c("Estimate", "Std. Error", "t.value", "p.value")
  tcv <- qt( 1-0.05/2, object$n-pp1)
coef <- cbind(`Coef`=est, `(95%`=est-tcv*se, `CI)`=est+tcv*se, `SE`=se, `p`=pval)
rownames(coef) <- object$xnames

ans <- list(coefficients=coef,fit=object)

class(ans) <- list('summary.bigRfit')

ans

}
