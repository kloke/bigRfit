gettauBDF0<-function(scoreMat,ehat,n,p=0,delta=0.8,hparm=2,...) {

sciint<-as.double(mad(ehat))

scrs <- rep(as.double(scoreMat[,scrs]),times=scoreMat[,count])
scrsD <- rep(as.double(scoreMat[,scrsD]),times=scoreMat[,count])

	.Fortran('nscale',
		as.integer(n),
		as.double(.Machine$double.eps^0.25),
		as.double(delta),
		as.double(sciint),
		as.integer(0),
		as.integer(p),
		as.double(ehat),
		as.integer(order(ehat)),
		as.double(scrs),
		as.double(rep.int(0,n)),
		as.double(rep.int(0,n)),
		tauhat=as.double(0),
		as.double(rep.int(0,5)),
		as.integer(0),
		as.integer(1000),
		as.double(0),
		as.double(scrsD),
		as.double(hparm),
		PACKAGE='Rfit'
	)$tauhat
		
		
}

