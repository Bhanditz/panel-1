dclike.panel<-function(indata, ainv, amat, evalues, ncov, nstage)
{
#Copyright 1994 Robert Gentleman
	nsubj <- length(indata)
	like <- 0
	for(ct in 1:nsubj) {
		tdata <- indata[[ct]]
		storage.mode(tdata$stage) <- "integer"
		storage.mode(tdata$time) <- "double"
		storage.mode(tdata$cov) <- "integer"
		if(tdata$len > 1) {
			likefun <- .Fortran("complike",
				as.vector(tdata$time),
				as.vector(tdata$stage),
				as.vector(tdata$cov),
				like = as.double(like),
				as.array(amat),
				as.array(ainv),
				as.matrix(evalues),
				as.integer(nstage),
				as.integer(ncov),
				as.integer(tdata$len))
			like <- like + likefun$like
		}
	}
	return(like)
}
