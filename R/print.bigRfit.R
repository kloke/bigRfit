print.bigRfit <- function(x,...) { 

  cat("Call:\n")
  print(x$call)
  cat("\nSample size = ",x$n,"\n")
  cat("Number of bins = ",x$B,"\n")

  cat("\nNumber of iterations:",x$iter, "\n")
  cat("Converge:",x$converge, "\n")

  cat("Estimates:\n")
  if(is.null(x$coef)) {
    if(!is.null(x$betahat0)) {
#      cat("Regression Coefficients:\n")
      cat("  Intercept (beta0) :\n")
      print(x$betahat0,...)
    } else {
#      cat("Regression Coefficients (Centered Design Matrix):\n")
       if(!is.null(x$alphahat)) {
         cat("  Intercept (alphahat) :\n")
          print(x$alphahat,...)
       }
    }
    cat("betahat:\n")
    print(x$betahat,...)
  } else {
    cat("\nCoefficients:\n")
    print(x$coef,...)
  }

}
