\name{annotTrack}
\alias{annotTrack}
\docType{data}
\title{Example Annotation Track}
\description{
  Example annotation track.
}
\usage{data(annotTrack)}

\format{
  A data frame with 1001 rows, each representing one annotation and 9 columns from standard gtf format. 
}

\details{
  This is an example annotation track as it was imported from the Ensemble ftp download page. In total there are 1001 annotations from the human genome in release 70.
}

\source{
  http://www.ensembl.org/info/data/ftp/index.html
}

\examples{
data(annotTrack)
head(annotTrack)
}
\keyword{datasets}
 
