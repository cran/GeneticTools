\name{plot.eqtl}
\docType{methods}
\alias{plot.eqtl}
\alias{plot-method}
\alias{plot,eqtl-method}

\title{Plot an eqtl Object}
\description{
 The function offers informative plots for an \code{eqtl} object.
}
\usage{
 \method{plot}{eqtl}(x, file=NULL, which=NULL, sig=0.01, ...)
}

\arguments{
  \item{x}{Object of class \code{eqtl}.}
  \item{file}{Store set of graphics under that file name.}
  \item{which}{Specifies for which genes should the plot be created.}
  \item{sig}{Chosen significance level.}
  \item{...}{Additional plotting parameters.}
}

\details{
 This function plots an \code{eqtl} object and all surrounding SNPs.
}

\author{Daniel Fischer}

\examples{
# Perform eQTL (single location, one gene):
# myeqtl <- eQTL(geneMatrix,genoData,locs,genoSamples,singleSamples,windowSize,method="LM")

# plot(myeqtl)

}

\keyword{methods}
\keyword{hplot}
