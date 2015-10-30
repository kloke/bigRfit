cppFunction('double tauhatCpp(NumericVector e, double h) {
  int n = e.size();
  double est = 0.0;
  double sqrt12 = sqrt(12);
  double c1 = sqrt12/h;

  for( int i=0; i < n; ++i ) {
    for( int j=0; j < n; ++j ) {
      est += c1 * R::dnorm( (e[i] - e[j])/h,0.0,1.0,0 );
    }
  }

  est = est/n/n;
  est = 1/est;
  return est;

}')


tauhat <- function(e,bw=bw.nrd(e)) {
  tauhatCpp(e,bw)
}


