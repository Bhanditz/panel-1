eddcmp<-function(inmat)
{
# Copyright 1994 Robert Gentleman
	if(!is.matrix(inmat))
		stop("eddcmp can only be used for matrices")
	n <- nrow(inmat)
	if(n!=ncol(inmat))
		stop("eddcmp can only be used for square matrices")
	emat <- matrix(0, nrow = n, ncol = n)
	eval1 <- vector(mode = "double", length = n)
	eval2 <- vector(mode = "double", length = n)
	tivec <- vector(mode = "integer", length = n)
	trvec <- vector(mode = "double", length = n)
	ierr <- 0
	z <- .Fortran("eddcmp",
		as.matrix(inmat),
		as.matrix(emat),
		as.vector(eval1),
		as.vector(eval2),
		as.integer(n),
		as.vector(tivec),
		as.vector(trvec),
		as.integer(ierr))
	if(z[[8]]!=0)
		warning("maximum number of iterations exceeded in fortran")
	if(max(abs(z[[4]]))!=0)
		warning("there are complex eigenvalues")
	if(sum(duplicated(z[[3]])) > 0)
		warning("some eigenvalues are duplicated")
	return(list(evectors = z[[2]], evalues = z[[3]], im.evalues = z[[4]]))
	
}
