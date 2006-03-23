\name{oligodata}
\alias{oligodata}
\alias{datA.rma}
\alias{datB.rma}
\alias{datA.mas5}
\alias{datB.mas5}
\alias{datA.amp}
\alias{datB.amp}
\docType{data}
\title{Example data for package maCorrSample}
\description{
Example expression data to demonstrate the functionality of the package: two data sources A and B with 30 patients and 1000 genes each, for each of which we have RMA expression values, (logarithmized) MAS5 expression values, and MAS5 absent/present calls. 

Correspondingly, we have six data matrices whose name are constructed as \code{dat[A|B].[rma|mas5|amp]}.
}
\usage{data(oligodata)}
\format{
All matrices have genes as rows and samples as columns.
}
\source{
These are small anonymized excerpts from a real breast cancer data set. 
}
\seealso{\code{\link{CorrSample}}, \code{\link{plot.corr.sample}}}
\examples{
data(oligodata)
str(datA.rma)
str(datB.rma)
str(datA.mas5)
str(datB.mas5)
str(datA.amp)
str(datB.amp)
}
\keyword{datasets}
