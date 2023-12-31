\name{plot.corr.sample}
\alias{plot.corr.sample}
\alias{panel.corr.sample}
\title{Plot correlation of random pairs of genes}
\description{
\code{plot.corr.sample} provides the main functionality of package \code{maCorrPlot}: it plots the correlation of random pairs of genes against their variability. Systematic deviations of the plot from a constant zero indicate lack of normalization of the underlying expression matrix.

Formally, \code{plot.corr.sample} is the plotting method for objects of class \code{corr.sample} generated by \code{\link{CorrSample}}.

\code{panel.corr.sample} is the panel function that does the actual plotting work.
}
\usage{
plot.corr.sample(x, ..., cond, groups, grid = TRUE, refline = TRUE, xlog = TRUE,
                 scatter = FALSE, curve = FALSE, ci = TRUE, nint = 10,
				 alpha=0.95, length = 0.1, xlab="Standard Deviation")

panel.corr.sample(x, y, grid = TRUE, refline = TRUE, xlog = TRUE, 
                  scatter = FALSE, curve = FALSE, ci = TRUE, nint = 10, 
				  alpha=0.95, length = 0.1, col.line, col.symbol, ...)				 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y}{for \code{plot.corr.sample}, \code{x} is an object of class \code{corr.sample}, generated by function \code{CorrSample} that contains the pre-computed correlations and standard deviations for the random pairs of genes; for \code{panel.corr.sample}, \code{x} and \code{y} are the x- and y-components (or standard deviation and correlation) of the pairs of genes to be plotted in a specific panel.}
  \item{\dots}{either more objects of class \code{corr.sample} or plotting arguments passed to the underlying \code{xyplot}.}
  \item{cond}{either a vector or a list of vectors describing multiple objects of class \code{corr.sample}; ignored if only one such object (\code{x}) is specified. See Details and Examples.}
  \item{groups}{a vector or a list of vectors giving group membership for the random pairs of genes in the \code{corr.sample} objects to be plotted, resulting in multiple overlayed plots for each object. See Details and Examples.}
  \item{grid}{logical value indicating whether to draw a reference grid}
  \item{refline}{logical value indicaitng whether to draw a horizontal reference line a zero.}
  \item{xlog}{logical value indicating whether to use log-scale on the horizontal axis.}
  \item{scatter}{logical value indicaitng whether the plot the individual pairwise correlations.}
  \item{curve}{logical value indicating whether to fit a simple model for lack of fit to the correlations.}
  \item{ci}{logical value indicating whether to add confidence intervals.}
  \item{nint}{number of intervals into which to divide the horizontal axis for calculating average correlations.}
  \item{alpha}{the level of confidence to be plotted.}
  \item{length}{the length of the horizontal ticks indicating the ends of the confidence intervals (in inches).}
  \item{xlab}{the label for the horizontal axis.}
  \item{col.line, col.symbol}{graphical parameters that control the color of the correlation lines and the scatter plotting symbols}
}
\details{
The underlying plotting engine is \code{xyplot}, using \code{panel.corr.sample} as panel function, which also interprets most of the graphical parameters. Note that two kinds of arguments can be specified via \code{...}: First, an unlimited number of extra \code{corr.sample} objects, in case we want to display different expression measures for the same expression matrix, or compare different expression matrices, or both; this is somewhat similar to the behaviour of \code{boxplot.default}. Second, everything that does not inherit from \code{corr.sample} is passed on to \code{xyplot}, so in theory, the full range of lattice control options is available, as long as they do not conflixt with named arguments to \code{plot.corr.sample}, like \code{xlog} or \code{xlab}. 

Two mechanisms for comparisons within the same plot are available: First, as mentioned above, multiple \code{corr.sample} objects can be shown in the same graph, each within its own panel. If no \code{cond} is specified, these panels are just numbered in the order in which the objects appear in the arguments. Alternatively, one or two factors can be associated with each factor: in the first case, \code{cond} is just a vector with as many entries as \code{corr.sample} objects in the argument list; these entries are used to label the panels of the corresponding \code{corr.sample} objects. In the second case, \code{cond} is a list with two such vectors, and the objects are cross-classified according to both categories, and the panels are arranged in a row-column pattern reflecting this cross-classification, see Examples.

The other mechanism for graphical comparisons within the same plot is via \code{groups}, which draws different  correlation curves for different sub-groups of pairs of genes; the standard example is to classify pairs of genes according to their common or average score in regard to a quality control measure like the MAS5 presence calls, see Examples. These sub-groups are specified via \code{groups}; if there is only one \code{corr.sample} object in the function call (\code{x}), \code{groups} is just a vector with as many entries as there are random paris of genes in \code{x}. If several objects of class \code{corr.sample} have been specified in the function call, \code{groups} is a list of as many vectors as objects, where each vector has as many entries as the corresponding object has pairs of genes. 
}
\value{
A plot created by \code{xyplot}.
}
\references{
Ploner A, Miller LD, Hall P, Bergh J, Pawitan Y. Correlation test to assess low-level processing of high-density oligonucleotide microarray data. BMC Bioinformatics, 2005, 6(1):80 
\url{http://www.pubmedcentral.gov/articlerender.fcgi?tool=pubmed&pubmedid=15799785}
}
\author{Alexander Ploner \email{Alexander.Ploner@ki.se}}
\section{Warning }{\code{cond} is translated into conditioning variables for \code{xyplot}, which will not hesitate to average correlations across different \code{corr.sample} objects. It's hard to see when this would be a good idea, therefore \code{plot.corr.sample} will generate a warning.}
\seealso{\code{\link{CorrSample}}, \code{\link[lattice]{xyplot}}}
\examples{
# Get small example data
data(oligodata)
dim(datA.rma)
dim(datB.rma)

# Compute the correlations for 500 random pairs, 
# Larger numbers are reasonable for larger data sets
cs1.rma = CorrSample(datA.rma, 500, seed=210)
plot(cs1.rma)

# Change the plot
plot(cs1.rma, scatter=TRUE, curve=TRUE, alpha=0.99)

# Compare with MAS5 values for the same data set
cs1.mas5 = CorrSample(datA.mas5, 500, seed=210)
plot(cs1.rma, cs1.mas5, cond=c("RMA","MAS5"))

# We group pairs of gene by their average number of MAS5 present calls
pcntA = rowSums(datA.amp[cs1.mas5$ndx1, ]=="P") +
        rowSums(datA.amp[cs1.mas5$ndx2, ]=="P")
hist(pcntA)
pgrpA = cut(pcntA, c(0, 20, 40, 60), include.lowest=TRUE)
table(pgrpA)

# Plot the RMA values according to their MAS5 status 
# The artificial correlation is due to gene pairs with few present calls
plot(cs1.rma, groups=pgrpA, nint=5, auto.key=TRUE, ylim=c(-0.3, 0.5))

# Combine grouping and multiple conditions
plot(cs1.rma, cs1.mas5, cond=c("RMA","MAS5"), groups=list(pgrpA, pgrpA), 
     nint=5, auto.key=TRUE, ylim=c(-0.3, 0.5))

# Compare with second data set
# Specify more than one condition
cs2.rma  = CorrSample(datB.rma, 500, seed=391)
cs2.mas5 = CorrSample(datB.mas5, 500, seed=391)
plot(cs1.rma, cs1.mas5, cs2.rma, cs2.mas5,     
     cond=list(c("RMA","MAS5","RMA","MAS5"), c("A","A","B","B")))

}

\keyword{hplot}% at least one, from doc/KEYWORDS

