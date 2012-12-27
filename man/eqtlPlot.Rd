\name{eqtlPlot}
\alias{eqtlPlot}

\title{
  Plot Expressions for Genotype groups
}

\description{
  This function plots expression value boxplots, grouped for the different genotypes.
}
\usage{
  eqtlPlot(snp=NULL, gene=NULL, eqtl=NULL,geno=NULL,gex=NULL, genoSamples=NULL)
}

\arguments{
  \item{snp}{Grouping SNP.}
  \item{gene}{Expression values of a gene.}
  \item{eqtl}{An eQTL object.}
  \item{geno}{An genotype object.}
  \item{gex}{An gene expression object.}
  \item{genoSamples}{Sample names of the genotype object.}
}

\details{
 This function can be used to plot the expression values of a certain gene and grouped with one SNP.
 
 This can be done after an eQTL is performed. In that case just the \code{eqtl} object has to be given. If
 no eQTL has been performed previous than the plot can be created be giving the objects \code{geno}, \code{gex}
 and \code{genoSamples} instead, see \code{\link{eQTL}}.

 The \code{gene} options gives the name of the gene whos expression values should be used and \code{snp} whos 
 SNPs genotypes.

 Currently only the eQTL option is implemented.
}
\value{
 A Figure with three boxplots, giving the expressions for the three genotype groups.
}

\seealso{
  \code{\link{eQTL}}
}

\author{
  Daniel Fischer
}

\examples{
# First perform an eQTL
# lm.myEQTL <- eQTL(gex=ourExpression, geno=genotData, xAnnot=xAnnotDF, windowSize=1, genoSamples=genoSamples, method="LM")

# eqtlPlot(snp="RS1928123", gene="ARLTS", eqtl=lm.myEQTL)
}

\keyword{
  graphic
}
 
