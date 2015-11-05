cppFunction('double tauhatCppGS(NumericVector e, NumericVector s, NumericVector c, double h, int n) {
  int m = e.size();
  double est = 0.0;

  for( int i=0; i < m; ++i ) {
    for( int j=0; j < m; ++j ) {
      est += s[i] * c[i]*c[j]*R::dnorm( (e[i] - e[j])/h,0.0,1.0,0 )/h;
    }
  }

  est = est/n/n;
  est = 1/est;
  return est;

}')


tauhat.gs <- function(mids,counts,scores,bw=bw.nrd(mids)) {
  tauhatCppGS(sort(mids),getScoresDeriv(scores,seq_len(length(mids))/length(mids)),counts,bw,sum(counts))
}


