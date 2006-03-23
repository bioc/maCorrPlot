\name{CorrSample}
\alias{CorrSample}
\alias{RandPairs}
\title{Sample correlations for random pairs of genes}
\description{
\code{CorrSample} calculates the correlations, standard deviations and some auxiliary variables for random pairs of genes. A plot of the resulting object that shows that these correlations dependend systematically on the genes' variability, suggests a lack of normalization.

\code{RandPairs} is a helper function for generating random pairs from a list of genes.
}
\usage{
CorrSample(x, np, seed, rp, ndx)

RandPairs(probes, number)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a gene expression matrix, with samples as columns and genes as rows; missing values are accepted.}
  \item{np, number}{the number of random pairs}
  \item{seed}{an optional seed for the random sampling}
  \item{rp}{an optional matrix with two columns specifying the random pairs, see Details.}
  \item{ndx}{an optional logical matrix of the same dimension as \code{x} that allows to eliminate a subset of the expression values from the calculation of the correlations, standard deviations and auxiliary variables.}
  \item{probes}{a vector of genes from which to draw random pairs; can be integer, as a vector of row indices, or character, as a vector of row names.}
}
\details{
The sample of random pairs can be specified in a replicable manner either via \code{np} and \code{seed}, or by using the output from \code{RandPairs} for the parameter \code{rp}. In case we want to use the same set of random pairs (e.g. when comparing different expression measures on the same data set), the second option will be faster.
}
\value{
An object of class \code{corr.sample}; this is just a data frame with an extra class tag to allow for a plotting method. 

The data frame has \code{np} rows and nine columns: 
\itemize{
\item{\code{Correlation}}{the correlation between the two genes across samples}
\item{\code{StdDev}}{the geometric mean of the standard deviations of the two genes}
\item{\code{sd1},\code{sd2}}{the standard deviations of the genes}
\item{\code{m1},\code{m2}}{the means of the genes}
\item{\code{ndx1},\code{ndx2}}{the indices of the two genes; by default, these will be the corresponding row indices of \code{x}, but if \code{rp} is specified, they might be gene names.}
}}
\references{
Ploner A, Miller LD, Hall P, Bergh J, Pawitan Y. Correlation test to assess low-level processing of high-density oligonucleotide microarray data. BMC Bioinformatics, 2005, 6(1):80 
\url{http://www.pubmedcentral.gov/articlerender.fcgi?tool=pubmed&pubmedid=15799785}
}
\author{Alexander Ploner \email{Alexander.Ploner@ki.se}}
\seealso{\code{\link{plot.corr.sample}}}
\examples{
# Get small example data
data(oligodata)
dim(datA.rma)

# Compute the correlations for 500 random pairs, 
# that is ca. 1/1000 of all possible pairs 
# Larger numbers are reasonable for larger data sets
cs1 = CorrSample(datA.rma, 500, seed=210)
cs1[1:5,]

# Clear correlation for pairs of genes with low average variability 
plot(cs1)

# A different way of specifying the same 
set.seed(210)
rp = RandPairs(rownames(datA.rma), 500)
cs2 = CorrSample(datA.rma, rp=rp)
cs2[1:5,]
plot(cs2)
}
\keyword{datagen}% at least one, from doc/KEYWORDS

