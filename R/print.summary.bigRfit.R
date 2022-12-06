print.summary.bigRfit <- function(x,digits=getOption("digits")-3,...) { 

  cat("Call:\n")
  print(x$fit$call)
  cat("\nSample size = ",x$fit$n,"\n")
  cat("Number of bins = ",x$fit$B,"\n")

  cat("\nNumber of iterations:",x$fit$iter, "\n")
  cat("Converge:",x$fit$converge, "\n")

  cat("Estimates:\n")
  print(round(x$coef,digits=digits),...)

}
