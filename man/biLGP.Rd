\name{biLGP}
\alias{biLGP}
\alias{dbiLGP}
\alias{biLGP.logMV}
\alias{rbiLGP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The bivariate Lagrangian Poisson (LGP) distribution}

\description{
Density, random-number generation, and moments of the log-transformed distribution.
}
\usage{
dbiLGP(y, theta, lambda, nc=NULL, log=FALSE, add.carefully=FALSE)
biLGP.logMV(theta,lambda,nc=NULL,const.add=1,tol=1e-14,add.carefully=FALSE)
rbiLGP(n, theta, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Numeric vector or two-column matrix of bivariate data.  If matrix, each row corresponds to an observation.}
  \item{theta}{Numeric vector or three-column matrix of non-negative values for index parameters \eqn{\theta _0}{theta0}, \eqn{\theta _1}{theta1}, and \eqn{\theta _2}{theta2}, in that order.  If matrix, is read by row.}
  \item{lambda}{Numeric vector or three-column matrix of values for multiplicative parameters \eqn{\lambda _0}{lambda0}, \eqn{\lambda _1}{lambda1}, and \eqn{\lambda _2}{lambda2}, in that order.  If matrix, is read by row.  Values must be on the interval [-1,1].}
  \item{nc}{Numeric vector or three-column matrix of (reciprocals of) the normalizing constants.  These constants differ from 1 only if the corresponding \code{lambda} parameter is negative; see \code{\link{dLGP}()} for details. If matrix, is read by row. Defaults to \code{NULL}, in which case the normalizing constants are computed automatically.}
  \item{log}{Logical; should the natural log of the probability be returned?  Defaults to \code{FALSE}.}
  \item{add.carefully}{Logical.  If \code{TRUE}, the program takes extra steps to try to prevent round-off error during the addition of probabilities.  Defaults to \code{FALSE}, which is recommended, since using \code{TRUE} is slower and rarely makes a noticeable difference in practice.}
  \item{const.add}{Numeric vector of positive constants to add to the non-negative integers before taking their natural logarithm.  Defaults to 1, for the typical \eqn{\log (y+1)}{log(y+1)} transformation.}
  \item{tol}{Numeric; must be positive.  When \code{biLGP.logMV()} is calculating the second moment of the log-transformed distribution, it stops when the next term in the series is smaller than \code{tol}.}
  \item{n}{Integer; number of observations to be randomly generated.}
}
\details{
The bivariate LGP is constructed from three independent latent random variables, \eqn{X_0}{X0}, \eqn{X_1}{X1}, and \eqn{X_2}{X2}, where
\deqn{X_0 \sim LGP(\theta _0, \lambda _0)}{X0 ~ LGP(theta0, lambda0)}
\deqn{X_1 \sim LGP(\theta _1, \lambda _1)}{X1 ~ LGP(theta1, lambda1)}
\deqn{X_2 \sim LGP(\theta _2, \lambda _2)}{X2 ~ LGP(theta2, lambda2)}
The observable variables, \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2}, are defined as \eqn{Y_1 = X_0 + X_1}{Y1 = X0 + X1} and \eqn{Y_2 = X_0 + X_2}{Y2 = X0 + X2}, and thus the dependence between \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2} arises because of the common term \eqn{X_0}{X0}.  The joint PMF of \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2} is derived from the joint PMF of the three independent latent variables, with \eqn{X_1}{X1} and \eqn{X_2}{X2} re-expressed as \eqn{Y_1 - X_0}{Y1 - X0} and \eqn{Y_2 - X_0}{Y2 - X0}, and after \eqn{X_0}{X0} is marginalized out.

Function \code{dbiLGP()} is the bivariate LGP density (PMF).  Function \code{rbiLGP()} generates random draws from the bivariate LGP distribution, via calls to \code{\link{rLGP}()}.  Function \code{biLGP.logMV()} numerically computes the means, variances, and covariance of a bivariate LGP distribution, after it has been log transformed following addition of a positive constant.

Vectors of numeric arguments other than \code{tol} are cycled, whereas only the first element of logical and integer arguments is used.
}
\value{
\code{dbiLGP()} returns a numeric vector of probabilities.  \code{rbiLGP()} returns a matrix of random draws, which is of type 'numeric' (rather than 'integer', even though the bivariate LGP only has support on the non-negative integers).  \code{biLGP.logMV()} returns a numeric matrix with the following five named columns:
  \enumerate{
  \item \code{EY1}: Post-tranformation expectation of \eqn{Y_1}{Y1}.
  \item \code{EY2}: Post-tranformation expectation of \eqn{Y_2}{Y2}.
  \item \code{VY1}: Post-tranformation variance of \eqn{Y_1}{Y1}.
  \item \code{VY2}: Post-tranformation variance of \eqn{Y_2}{Y2}.
  \item \code{COV}: Post-tranformation covariance of \eqn{Y_1}{Y1} and \eqn{Y_2}{Y2}.
  }
}

\references{
Famoye, F., & Consul, P. C.  (1995).  Bivariate generalized Poisson distribution with some applications.  \emph{Metrika}, 42, 127-138.

Consul, P. C., & Famoye, F.  (2006).  \emph{Lagrangian Probability Distributions}.  Boston: Birkhauser.
}
\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}
}

\seealso{
\code{\link{LGP}}, \code{\link[stats:Poisson]{dpois}()}
}
\examples{
## The following two lines do the same thing:
dbiLGP(y=1,theta=1,lambda=0.1)
dbiLGP(y=c(1,1),theta=c(1,1,1),lambda=c(0.1,0.1,0.1))

dbiLGP(y=c(1,1,2,2,3,5),theta=c(1,1,1,2,2,2),lambda=0.1)
## Due to argument cycling, the above line is doing the following three steps:
dbiLGP(y=c(1,1),theta=c(1,1,1),lambda=c(0.1,0.1,0.1))
dbiLGP(y=c(2,2),theta=c(2,2,2),lambda=c(0.1,0.1,0.1))
dbiLGP(y=c(3,5),theta=c(1,1,1),lambda=c(0.1,0.1,0.1))

## Inputs to dbiLGP() can be matrices, too:
dbiLGP(y=matrix(c(1,1,2,2,3,5),ncol=2,byrow=TRUE),
  theta=matrix(c(1,1,1,2,2,2,1,1,1),ncol=3,byrow=TRUE),
  lambda=0.1)

## theta0 = 0 implies independence:
a <- dbiLGP(y=c(1,3),theta=c(0,1,2),lambda=c(0.1,-0.1,0.5))
b <- dLGP(x=1,theta=1,lambda=-0.1) * dLGP(x=3,theta=2,lambda=0.5)
a-b #<--near zero.
## lambdas of zero yield the ordinary Poisson:
a <- dbiLGP(y=c(1,3), theta=c(0,1,2),lambda=0)
b <- dpois(x=1,lambda=1) * dpois(x=3,lambda=2) #<--LGP theta is pois lambda
a-b #<--near zero

( y <- rbiLGP(10,theta=c(1.1,0.87,5.5),lambda=c(0.87,0.89,0.90)) )
dbiLGP(y=y,theta=c(1.1,0.87,5.5),lambda=c(0.87,0.89,0.90))

biLGP.logMV(theta=c(0.65,0.35,0.35),lambda=0.7,tol=1e-8)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
