\name{eQTL}
\alias{eQTL}

\title{
  Performs an eQTL Analysis.
}

\description{
  This function performs an eQTL analysis.
}
\usage{
  eQTL(gex, geno, xAnnot, xSamples=NULL, genoSamples, windowSize=0.5, method="LM")
}

\arguments{
  \item{gex}{Matrix/Vector with expression values.}
  \item{geno}{Genotype data.}
  \item{xAnnot}{Location annotations for the expression values.}
  \item{xSamples}{Sample names for the expression values, see details.}
  \item{genoSamples}{Sample names for the genotype values, see details.}
  \item{windowSize}{Size of the window around the center gene.}
  \item{method}{Method of choice for the eQTL, see details.}
}

\details{
This function performs an eQTL analysis with different types of tests. The type of test
can be specified with the \code{method} option. Possible options are \code{"LM"} and \code{"directional"}. 
\code{"LM"} fits for each SNP within a predefined window of size \code{windowSize} (in MB) around a gene
a linear model for the genotype information and the corresponding gene expression. The tested hypothesis
is then if the slope is equal to zero or not.

The \code{"directional"} option applies a directional test based on probabilistic indices for triples as described 
in Fischer et al (2012). The test is applied for the two probabilistic indices \eqn{P_{0,1,2}} and
\eqn{P_{2,1,0}} and we combine the two corresponding p-values \eqn{p_{012}=p_1} and \eqn{p_{210}=p_2} from previous tests then as overall p-value
\eqn{min(2 min(p_1 , p_2 ), 1)}. We refer here to the different genotype groups as 0,1,2.

The gene expressions are given in \code{gex}. If several genes should be tested \code{gex} is a matrix and each
column referes to a gene and each row to an individuum. The column names of this matrix should match then with the 
names used in the \code{xAnnot} object. Sample names can either be given as rownames in the matrix or as separate
vector in \code{xSamples}. If only one gene should be tested then \code{gex} can be a vector.

The genotype information is provided in the \code{geno} object. Here one can either specify the file
name of a ped/map file pair. In that case the function imports the genotype information using the
\code{SnpStats} package. In case the genotype information has been imported already earlier using
\code{SnpStats::read.pedfile()} the resulting \code{SnpMatrix} can also given as a parameter for \code{geno}.

The \code{xAnnot} object carries the annotations for the gene expressions. In case of multiple locations per gene 
it is of type list and each list item stores the information for one gene. Within the list item is then a data frame
with the three columns \code{Chr}, \code{Start}, \code{End} and each row refers to one matching chromosomal postion of the underlying gene.
Especially when probes of ssRNA are considered the chromosomal positions of a probe are not necessary unique. The names
of the list \code{xAnnot} are the names of the genes and they have to match with the column names of \code{gex}. However, the order
does not have to be the same, and \code{xAnnot} can include more genes than given in \code{gex}. The function finds
and uses then the union between the column names of \code{gex} and the list entries of \code{xAnnot}.
Alternative \code{xAnnot} can be a data frame if no multiple locations are considered. In that case \code{xAnnot} has
to be a data frame with the four columns \code{SNP}, \code{Chr}, \code{Start}, \code{End}.

The option \code{genoSamples} is used in case that the sample names in the ped/map file (or SnpMatrix)
do not match with \code{rownames(gex)} given in the expression matrix. The vector \code{genoSamples} is as
long as the \code{geno} object has samples, but gives then for each individuum in geno the corresponding
name in the \code{gex} object. The function finds then also the smallest union between the two data objects.
If there are repeated measurements per individual for the genotypes we take by default only the first
appearance in the data and neglect all successive values. Currently this cannot be changed. In case
this behavior is not desired, the user has to remove the corresponding rows from \code{geno}.

}
\value{
 A list of class \code{eqtl} containing the values 
 \item{gex}{The \code{gex} object from the function call.}
 \item{geno}{The \code{geno} object from the function call.}
 \item{xAnnot}{The \code{xAnnot} object from the function call.}
 \item{genoSamples}{The \code{genoSamples} object from the function call.}
 \item{windowSize}{The \code{windowSize} object from the function call.}

 and an incapsulated list \code{eqtl} where each list item is a tested gene location and contains the items

 \item{ProbeLoc}{Used position of that gene. (Only different from 1 if multiple locations are considered.)}
 \item{TestedSNP}{Details about the considered SNPs.}
 \item{p.values}{P values of the test.}
 \item{GeneInfo}{Details about the center gene.}
}

\references{
  \cite{Fischer, D., Oja, H., Sen, P.K., Schleutker, J., Wahlfors, T. (2012): Generalized Mann-Whitney Type Tests for Microarray Experiments, submitted article.}

  \cite{Fischer, D., Oja, H. (2012): A permutation type test for calculating generalized Mann-Whitney tests, submitted article.}
}


\author{
  Daniel Fischer
}

\examples{
# Will be filled as soon as data is available (the corresponding article is currently
# under review). The data will then be published as a accompanying article.

# Import first the expression data
# geneMatrix <- 
# geneVector <- 

# Then import the genotype information in map/ped format
# genoData <- 

# Load the expression annotations
# multLoc <- 
# singleLoc <- 

# Optional the sample names
# genoSamples <- 
# expSamples <- 

# Set the window size in MB (here then 500.000 base pairs to BOTH sides:)
# windowSize <- 0.5

# Perform different eQTLS:

# eqtl1 <- eQTL(geneMatrix, genoData, multLoc, genoSamples, expSamples, windowSize)
# eqtl2 <- eQTL(geneMatrix, genoData, singleLoc, genoSamples, expSamples, windowSize)

}

\keyword{
  graphic
}
 
