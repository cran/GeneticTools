\name{calcDistances}
\alias{calcDistances}

\title{
  Calculate Distances Between Time Series
}

\description{
  This is an experimental function to calculate disances between time series and shouldn't be used in productive code.
}
\usage{
  calcDistances(X, method="C", nodes=20)
}

\arguments{
  \item{X}{Matrix with time series values.}
  \item{method}{Method  of choice to calculate the distance.}
  \item{nodes}{Number of splines for the spline-based approach.}
}

\details{
 This is an early attempt to implement a distance measure for time series clustering. It has not yet been tested
 and hence should be used carefully.

 Currently there are two different methods available, the default is \code{method="C"}. The alternative is \code{spline},
 what is a spline based approach.

 The matrix \code{X} contains the individuals in the rows and the columns represent the time points. 
 
 As this is experimental code, we do not give further information, yet.  
}
\value{
 A matrix with pairwise distances.
}

\author{
  Daniel Fischer
}

\keyword{
  methods
}
 
