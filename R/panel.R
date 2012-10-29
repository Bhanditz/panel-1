panel<- function(indata, qmatf, gamma, qderivf, npar, nstage, ncov, 
verbose = FALSE, tol = 0.001)
{
#Copyright 1994 Robert Gentleman
	stepsize <- 1
	amat <- array(0, c(ncov, nstage, nstage))
	ainv <- array(0, c(ncov, nstage, nstage))
	evalues <- matrix(0, nrow = ncov, ncol = nstage)
	garray <- array(0, dim = c(ncov, npar, nstage, nstage))
	for(iter in 1:40) {
		wqm <- qmatf(gamma)
		wderivs <- qderivf(gamma)
		for(i in 1:ncov) {
			qdcmp <- eddcmp(wqm[i,  ,  ])
			amat[i,  ,  ] <- qdcmp$evectors
			ainv[i,  ,  ] <- solve(qdcmp$evectors)
			evalues[i,  ] <- qdcmp$evalues
			for(j in 1:npar)
				garray[i, j,  ,  ] <- ainv[i,  ,  ] %*% wderivs[
				  j, i,  ,  ] %*% amat[i,  ,  ]
		}
		nsubj <- length(indata)
		info <- matrix(0, nrow = npar, ncol = npar)
		score <- rep(0, npar)
		for(tdata in indata) {
			if(tdata$len > 1) {
				storage.mode(tdata$stage) <- "integer"
				storage.mode(tdata$time) <- "double"
				storage.mode(tdata$cov) <- "integer"
				tcmp <- .Fortran("cmpscore",
				  as.vector(tdata$time),
				  as.vector(tdata$stage),
				  as.vector(tdata$cov),
				  as.array(garray),
				  as.array(amat),
				  as.array(ainv),
				  as.matrix(evalues),
				  as.integer(tdata$len),
				  score = as.vector(score),
				  info = as.matrix(info),
				  as.vector(gamma),
				  as.integer(ncov),
				  as.integer(npar),
				  as.integer(nstage),
				  PACKAGE = "panel")
				score <- tcmp$score
				info <- tcmp$info
			}
		}
		if(verbose) {
			cat("Iteration")
			print(iter)
			cat("gamma", fill = TRUE)
			print(gamma)
			cat("score at end of iteration", fill = TRUE)
			print(score)
		}
		if(max(abs(score)) < tol) break	#no point in going further
		rightstep <- FALSE
		dir <- solve(info, score)
		lval1 <- dclike.panel(indata, ainv, amat, evalues, ncov, nstage
			)
		minstep <- 0.0001
		while((!(rightstep)) && (stepsize > minstep)) {
			nstep <- stepsize * dir
			ngamma<- gamma + nstep
			pred <- score %*% nstep - nstep %*% info %*% nstep/2
			q2arr <- qmatf(ngamma)
			for(i in 1:ncov) {
				qdcmp <- eddcmp(q2arr[i,  ,  ])
				amat[i,  ,  ] <- qdcmp$evectors
				ainv[i,  ,  ] <- solve(qdcmp$evectors)
				evalues[i,  ] <- qdcmp$evalues
			}
			lval2 <- dclike.panel(indata, ainv, amat, evalues, ncov,
				nstage)
			obsvd <- lval2 - lval1
			if(obsvd/pred < 0.25)
				stepsize <- stepsize/2
			else if(obsvd/pred > 0.75) {
				rightstep <- TRUE
				stepsize <- min(stepsize * 2, 1)
			}
			else rightstep <- TRUE
		}
		if(verbose) {
			ss <- paste("Log Likelihood ", lval2, " stepsize ", 
				stepsize)
			print(ss)
		}
		if(rightstep)
			gamma <- ngamma
		else stop("panel: no step in search direction possible")
	}
	cat(" ", fill = TRUE)
	cat("Results at convergence", fill = TRUE)
	cat("----------------------", fill = TRUE)
	cat("Log Likelihood")
	print(lval2)
	cat("gamma", fill = TRUE)
	print(gamma)
	cat("Score", fill = TRUE)
	print(score)
	cat("stderr(gamma)", fill = TRUE)
	print(sqrt(diag(solve(info))))
	return(list(gamma = gamma, info = info, like = lval2))
}
