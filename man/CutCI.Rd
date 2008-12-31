\name{CutCI}
\alias{CutCI}
\alias{CIrho}
\title{Calculate confidence intervals for grouped values}
\description{
\code{CutCI} groups values of one variable into intervals with the same number of observations each and computes confidence intervals for the mean of another variable in each interval.

\code{CIrho} computes the normal theory confidence interval for a vector of values.
}
\usage{
CutCI(dat, number = 10, func = mean, alpha=0.95)

CIrho(rho, alpha = 0.95)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{a numerical data frame or matrix with two columns, the first of which gets averaged, and the second of which defines the grouping}
  \item{number}{the number of equal-count intervals}
  \item{func}{summary function for computing the mean}
  \item{rho}{a vector of measurements}
  \item{alpha}{the desired confidence level}  
}
\details{
The quantiles for the confidence interval are taken from the standard normal distribution, so a reasonable number of observations per interval would be good. 
}
\value{
\code{CutCI} returns invisibly a list of length three:
  \item{x}{the midpoints of the grouping intervals}
  \item{y}{the means within each interval, as computed by \code{func}}
  \item{yci}{a matrix with two columns, giving the lower and upper end of the confidence interval respectively}


\code{CIrho} returns a vector of length two, containing the lower and upper end of the confidence interval.
}
\seealso{\code{\link{co.intervals}}}
\examples{
x = rnorm(100, mean=2)
CIrho(x)

y = 2 + 3*x + rnorm(100)
cc = CutCI(cbind(x,y), number=5)
print(cc)

# Show it
plot(cc$x, cc$y)
arrows(cc$x, cc$yci[,1], cc$x, cc$yci[,2], length=0)
}
\keyword{utilities}% at least one, from doc/KEYWORDS

