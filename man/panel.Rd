\name{panel}
\alias{panel}
\title{
A function to compute the maximum likelihood estimates of the transition
parameters from panel data.
}
\description{
The algorithm is discussed in Kalbfleisch and Lawless (1985) and
further in Gentleman (1994).
}
\usage{
panel(indata, qmatf, gamma, qderivf, npar, nstage, ncov, verbose = F,
        tol = 0.001)
}
\arguments{
\item{indata}{
indata is a list with one element for each individual. Each element of
indata is a structure with 4 components. 
time: a vector containing the inspection times
stage: a vector containing the state of the individual at the 
corresponding inspection time 
cov: a vector containing the covariate values at the inspection times
len: the number of inspections on 
that individual
}
\item{npar}{
npar is the number of parameters that are to be estimated
}
\item{nstage}{
nstage is the number of stages/states in the model
}
\item{ncov}{
ncov is the number of levels the covariate comes at
}
\item{qmatf}{
qmatf is a function that accepts as input gamma (the parameters)
and returns an array (ncov x nstage x nstage) with the ncov transition
intensity matrices.
}
\item{gamma}{
gamma is a vector (npar long) of estimates of the transition parameters.
}
\item{qderivf}{
qderivf is a function that takes gamma as an argument and returns an
array (npar x ncov x nstage x nstage) with the derivatives of the transition
intensity matrices with respect to each of the parameters.
}
\item{verbose}{
If verbose is set to TRUE then on each iteration the current estimate of theta
and the current value of the score function are printed out.
}}
\value{
panel returns a list with the parameter estimates, their standard errors
and the log likelihood on convergence.
}
\section{Side Effects}{
There are no known side effects.
}
\references{
Kalbfleisch J. D. and Lawless J. F. "The Analysis of Panel Data Under
a Markov Assumption", JASA (1985), p. 863.
Gentleman R. (1994). The Use of Covariate Information in Multi-state
Markov Models. Technical Report 13, Department of Statistics, University
of Auckland.
}
\seealso{
eddcmp, dclike.panel
}
\examples{
#See the accompanying document "A Users Guide to panel" for
#examples and code segments.
}
\keywords{model}







