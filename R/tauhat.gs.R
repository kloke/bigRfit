tauhat.gs <-
function(mids,counts,scores,bw=bw.nrd(mids)) {
  .Call('bigRfit_tauhatCppGS',mids,getScoresDeriv(scores,seq_len(length(mids))/(length(mids)+1)),counts,bw,sum(counts),PACKAGE='bigRfit')
}
