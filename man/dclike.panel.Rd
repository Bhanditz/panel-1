\name{dclike.panel}
\alias{dclike.panel}
\title{
A function to compute the value of the log likelihood for a given set of
panel data and a given value of the parameters.
}
\description{
A fortran program is called to compute the value of the log likelihood.
}
\usage{
dclike.panel(indata, ainv, amat, evalues, ncov, nstage)
}
\arguments{
\item{indata}{
indata is a list with one entry for each individual. Each entry is a structure
with the components:
time: the inspection times for that individual
stage: the stage the individual is in at each inspection times
cov: the value of the covariate for that individual at the inspection time
len: the number of inspections for that individual
}
\item{ainv}{
The inverse of the intensity matrix.
}
\item{amat}{ The intensity matrix.}
\item{evalues}{The eigenvalues.}
\item{ncov}{
The number of levels the covariate comes at.
}
\item{nstage}{
nstage is the number of stages in the Markov model.
}}
\value{
The value of the log likelihood is returned.
}
\section{Side Effects}{
There are no known side effects.
}
\seealso{
panel, eddcmp
}
\examples{
}
\keyword{utilities}

