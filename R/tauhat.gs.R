tauhat.gs <-
function(mids,counts,scores,bw=bw.nrd(mids)) {
  .Call('bigRfit_tauhatCppGS',sort(mids),getScoresDeriv(scores,seq_len(length(mids))/length(mids)),counts,bw,sum(counts),PACKAGE='bigRfit')
}
