\name{eddcmp}
\alias{eddcmp}
\title{Eigenvalue--eigenvector decomposition.}
\description{
A function that produces an eigenvalue eigenvector decomposition of a
transition intensity matrix for a Markov chain in double precision.
}
\usage{
eddcmp(inmat)
}
\arguments{
\item{inmat}{
inmat is any square matrix describing the markov chain
}}
\value{
A list with named components is returned. 

\item{evectors}{
 a matrix containing the left eigenvectors 
}
\item{evalues}{
a vector containing the eigenvalues 
}
\item{im.evalues}{
a vector which contains the imaginary components of the eigenvalues
if there are any (otherwise this is the zero vector).
}}
\examples{
}
\section{COMMENTS}{
  The matrix inmat is decomposed as evectors \%*\% diag(evalues) \%*\% solve(evectors).
For most matrices of interest the eigenvalues will all be real so the imaginary
part is commonly ignored.
The eigenvalues are ordered from smallest to largest.
}
\keyword{array}
