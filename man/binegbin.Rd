\name{binegbin}
\alias{dbinegbin}
\alias{binegbin.logMV}
\alias{rbinegbin}

\title{The bivariate negative binomial distribution}
\description{
Functions for the bivariate negative binomial distribution, as generated via trivariate reduction: density, random-number generation, and moments of the log-transformed distribution.
}
\usage{
dbinegbin(y, nu, p, log=FALSE, add.carefully=FALSE)
binegbin.logMV(nu,p,const.add=1,tol=1e-14,add.carefully=FALSE)
rbinegbin(n, nu, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Numeric vector or two-column matrix of bivariate data.  If matrix, each row corresponds to an observation.}
  \item{nu}{Numeric vector or three-column matrix of non-negative values for index parameters \eqn{\nu _0}{nu0}, \eqn{\nu _1}{nu}, and \eqn{\nu _2}{nu2}, in that order.  If matrix, is read by row.}
  \item{p}{Numeric vector or three-column matrix of values for Bernoulli parameters \eqn{p_0}{p0}, \eqn{p_1}{p1}, and \eqn{p_2}{p2}, in that order.  If matrix, is read by row.  Values must be on the interval (0,1].}
  \item{log}{Logical; should the natural log of the probability be returned?  Defaults to \code{FALSE}.}
  \item{add.carefully}{Logical.  If \code{TRUE}, the program takes extra steps to try to prevent round-off error during the addition of probabilities.  Defaults to \code{FALSE}, which is recommended, since using \code{TRUE} is slower and rarely makes a noticeable difference in practice.}
  \item{const.add}{Numeric vector of positive constants to add to the non-negative integers before taking their natural logarithm.  Defaults to 1, for the typical \eqn{\log (y+1)}{log(y+1)} transformation.}
  \item{tol}{Numeric; must be positive.  When \code{binegbin.logMV()} is calculating the second moment of the log-transformed distribution, it stops when the next term in the series is smaller than \code{tol}.}
  \item{n}{Integer; number of observations to be randomly generated.}
}
\details{
This bivariate negative binomial distribution is constructed from three independent latent variables, in the same manner as the \link[RMKdiscrete:biLGP]{bivariate Lagrangian Poisson} distribution.

Function \code{dbinegbin()} is the bivariate negative binomial density (PMF).  Function \code{rbinegbin()} generates random draws from the bivariate negative binomial distribution, via calls to \code{\link[stats:rnbinom]{rnbinom}()}.  Function \code{binegbin.logMV()} numerically computes the means, variances, and covariance of a bivariate LGP distribution, after it has been log transformed following addition of a positive constant.

Vectors of numeric arguments other than \code{tol} are cycled, whereas only the first element of logical and integer arguments is used.
}
\value{
\code{dbinegbin()} returns a numeric vector of probabilities.  \code{rbinegbin()} returns a matrix of random draws, which is of type 'numeric' (rather than 'integer', even though the negative binomial only has support on the non-negative integers).  \code{binegbin.logMV()} returns a numeric matrix with the following five named columns:
  \enumerate{
  \item \code{EY1}: Post-tranformation expectation of \eqn{Y_1}{Y1}.
  \item \code{EY2}: Post-tranformation expectation of \eqn{Y_2}{Y2}.
  \item \code{VY1}: Post-tranformation variance of \eqn{Y_1}{Y1}.
  \item \code{VY2}: Post-tranformation variance of \eqn{Y_2}{Y2}.
  \item \code{COV}: Post-tranformation covariance of \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2}.
  }
}

\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}
}

\seealso{
\code{\link{dbiLGP}}, \code{\link[stats:dnbinom]{dnbinom}()},  \code{\link[stats:rnbinom]{rnbinom}()}
}
\examples{
## The following two lines do the same thing:
dbinegbin(y=1,nu=1,p=0.9)
dbinegbin(y=c(1,1),nu=c(1,1,1),p=c(0.9,0.9,0.9))

dbinegbin(y=c(1,1,2,2,3,5),nu=c(1,1,1,2,2,2),p=0.9)
## Due to argument cycling, the above line is doing the following three steps:
dbinegbin(y=c(1,1),nu=c(1,1,1),p=c(0.9,0.9,0.9))
dbinegbin(y=c(2,2),nu=c(2,2,2),p=c(0.9,0.9,0.9))
dbinegbin(y=c(3,5),nu=c(1,1,1),p=c(0.9,0.9,0.9))

## Inputs to dbinegbin() can be matrices, too:
dbinegbin(y=matrix(c(1,1,2,2,3,5),ncol=2,byrow=TRUE),
  nu=matrix(c(1,1,1,2,2,2,1,1,1),ncol=3,byrow=TRUE),
  p=0.9)
  
## nu0 = 0 implies independence:
a <- dbinegbin(y=c(1,3),nu=c(0,1,2),p=c(0.1,0.5,0.9))
b <- dnegbin(x=1,nu=1,p=0.5) * dnegbin(x=3,nu=2,p=0.9)
a-b #<--near zero.

( y <- rbinegbin(10,nu=c(1.1,0.87,5.5),p=c(0.87,0.89,0.90)) )
dbinegbin(y=y,nu=c(1.1,0.87,5.5),p=c(0.87,0.89,0.90))
( mv <- negbinMVP(nu=c(1.1,0.87,5.5),p=c(0.87,0.89,0.90)) )
mv[1,2] #<--Covariance of this distribution
mv[1,2]+mv[2,2] #<--Marginal variance of Y1
mv[1,2]+mv[3,2] #<--Marginal variance of Y2
mv[1,2]/(sqrt(mv[1,2]+mv[2,2])*sqrt(mv[1,2]+mv[3,2])) #<--Correlation
logmv <- binegbin.logMV(nu=c(1.1,0.87,5.5),p=c(0.87,0.89,0.90))
## Log transformation nearly cuts the correlation in half:
logmv[1,5]/sqrt(logmv[1,3]*logmv[1,4])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
