\name{kldata}
\alias{kldata}
\non_function{}
\title{ Simulated panel data from "The analysis
of panel data under a Markov assumption"  }
\usage{data(kldata)}
\description{
   These are simulated data which consists of 120 observations.
  15 individuals were started in each of States 1 and 2 for each of four
groups with regression parameters $(z_1,z_2)=$ (0,0), (0,1), (1,0) and
(1,1).
Individuals were followed for 5 transitions at equidistant time
intervals.
}
\format{
The data are in the form of a list. Each element of the list consists of
four components.
\itemize{
  \item time 
  \item stage
  \item cov  
  \item len
  }
This would be a very good candidate for the \code{\link[lgtdl]{lgtdl}} class of
objects. 
}

\source{
Kalbfleisch, J.D. and Lawless, J.F. 'The analysis
of panel data under a Markov assumption', \emph{Journal of the 
American Statistical Association}, \bold{ 80}, 863-871, (1985).
}

\examples{
data(kldata)
}
\keyword{datasets}
