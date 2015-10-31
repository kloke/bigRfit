cppFunction('double tauhatCppGS(NumericVector e, NumericVector s, double h) {
  int n = e.size();
  double est = 0.0;

  for( int i=0; i < n; ++i ) {
    for( int j=0; j < n; ++j ) {
      est += s[i] * R::dnorm( (e[i] - e[j])/h,0.0,1.0,0 )/h;
    }
  }

  est = est/n/n;
  est = 1/est;
  return est;

}')


tauhat.gs <- function(e,scores,bw=bw.nrd(e)) {
  tauhatCppGS(sort(e),getScoresDeriv(scores,seq_len(length(e))/length(e)),bw)
}


