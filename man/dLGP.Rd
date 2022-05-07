\name{LGP}
\alias{LGP}
\alias{LGPMVP}
\alias{LGP.findmax}
\alias{LGP.get.nc}
\alias{dLGP}
\alias{pLGP}
\alias{qLGP}
\alias{rLGP}
\alias{sLGP}

\title{The (univariate) Lagrangian Poisson (LGP) Distribution}
\description{Density, distribution function, quantile function, summary, random number generation, and utility functions for the (univariate) Lagrangian Poisson distribution.}
\usage{
dLGP(x,theta,lambda,nc=NULL,log=FALSE)
pLGP(q,theta,lambda,nc=NULL,lower.tail=TRUE,log.p=FALSE,add.carefully=FALSE)
qLGP(p,theta,lambda,nc=NULL,lower.tail=TRUE,log.p=FALSE,add.carefully=FALSE)
rLGP(n,theta,lambda)
sLGP(theta,lambda,nc=NULL,do.numerically=FALSE,add.carefully=FALSE)
LGP.findmax(theta,lambda)
LGP.get.nc(theta,lambda,nctol=1e-14,add.carefully=FALSE)
LGPMVP(mu,sigma2,theta,lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, q}{Numeric vector of quantiles.}
  \item{p}{Numeric vector of probabilities.}
  \item{n}{Integer; number of observations to be randomly generated.}
  \item{theta}{Numeric; the index (or "additive") parameter of the LGP distribution.  Must be non-negative.}
  \item{lambda}{Numeric; the dispersion (or "Lagrangian" or "multiplicative") parameter of the LGP distribution.  Must not exceed 1 in absolute value.  When equal to zero, the LGP reduces to the ordinary Poisson distribution, with mean equal to \code{theta}.  When negative, then the distribution has an upper limit to its support, which may be found with \code{LGP.findmax()}}.
  \item{nc}{Numeric; the reciprocal of the normalizing constant of the distribution, by which the raw PMF must be multiplied so that it is a proper PMF, with values that sum to 1 across the support, when \code{lambda} is negative.  Defaults to \code{NULL}, in which case it is computed numerically by a call to \code{LGP.get.nc()}.}
  \item{log, log.p}{Logical; if \code{TRUE}, then probabilities p are given as log(p).}
  \item{lower.tail}{Logical; if \code{TRUE} (default), probabilities are \eqn{P[X \leq x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{nctol}{Numeric; while numerically computing the normalizing constant, how close to 1 should it be before stopping?  Ignored unless \code{lambda} is negative, and the upper support limit exceeds 200,000.}
  \item{add.carefully}{Logical.  If \code{TRUE}, the program takes extra steps to try to prevent round-off error during the addition of probabilities.  Defaults to \code{FALSE}, which is recommended, since using \code{TRUE} is slower and rarely makes a noticeable difference in practice.}
  \item{do.numerically}{Logical; should moments be computed numerically when \code{lambda}<0?  Defaults to \code{FALSE}, which is recommended unless the upper support limit is fairly small (say, less than 10).}
  \item{mu}{Numeric vector of mean parameters.}
  \item{sigma2}{"Sigma squared"--numeric vector of variance parameters.}
}
\details{
The Lagrangian Poisson (LGP) distribution has density
\deqn{p(x)=\frac{\theta (\theta + \lambda x)^{x-1} \exp(- \theta - \lambda x)}{x!}}{\theta*(\theta+\lambda*x)^(x-1)*exp(-\theta-\lambda*x)/x!} for \eqn{0,1,2,\ldots},
\deqn{p(x)=0} for \eqn{x>m} if \eqn{\lambda<0},
and zero otherwise, where \eqn{\theta>0}, \eqn{m=\lfloor-\theta / \lambda\rfloor} if \eqn{\lambda<0}, and \eqn{\max(-1,-\theta / m)\leq\lambda\leq 1}.  So, when \eqn{\lambda} is negative, there is an upper limit to the distribution's support, \eqn{m}, equal to \eqn{-\theta / \lambda}, rounded down to the next-smallest integer.  When \eqn{\lambda} is negative, the PMF must also be normalized numerically if it is to describe a proper probability distribution.  When \eqn{\lambda=0}, the Lagrangian Poisson reduces to the ordinary Poisson, with mean equal to \eqn{\theta}.  When \eqn{\theta=0}, we define the distribution as having unit mass on the event \eqn{X=0}.

Function \code{LGP.findmax()} calculates the value of upper support limit \eqn{m}; \code{LGP.get.nc()} calculates the (reciprocal of) the normalizing constant.

Function \code{LGPMVP()} accepts exactly two of its four arguments, and returns the corresponding values of the other two arguments.  For example, if given values for \code{theta} and \code{lambda}, it will return the corresponding means (\code{mu}) and variances (\code{sigma2}) of an LGP distribution with the given values of \eqn{\theta} and \eqn{\lambda}.  \code{LGPMVP()} does not enforce the parameter space as strictly as other functions, but will throw a warning for bad parameter values.

When the upper support limit is 5 or smaller, \code{rLGP()} uses simple inversion (i.e., random unit-uniform draws passed to \code{qLGP()}).  Otherwise, it uses random-number generation algorithms from Consul & Famoye (2006); exactly which algorithm is used depends upon the values of \code{theta} and \code{lambda}.  All four of \code{rLGP()}, \code{dLGP()}, \code{pLGP()}, and \code{qLGP()} make calls to the corresponding functions for the ordinary Poisson distribution (\code{\link[stats:dpois]{dpois}()},etc.) when \code{lambda}=0.

Vectors of numeric arguments are cycled, whereas only the first element of logical and integer arguments is used.

}
\value{
\code{dLGP()} and \code{pLGP()} return numeric vectors of probabilities.  \code{qLGP()}, \code{rLGP()}, and \code{LGP.findmax()} return  vectors of quantiles, which are of class 'numeric' rather than 'integer' for the sake of compatibility with very large values.  \code{LGP.get.nc()} returns a numeric vector of reciprocal normalizing constants.  \code{LGPMVP()} returns a numeric matrix with two columns, named for the missing arguments in the function call.

\code{sLGP()} returns a numeric matrix with 10 columns, with the mostly self-explanatory names \code{"Mean"}, \code{"Median"}, \code{"Mode"}, \code{"Variance"}, \code{"SD"}, \code{"ThirdCentralMoment"}, \code{"FourthCentralMoment"}, \cr \code{"PearsonsSkewness"}, \code{"Skewness"}, and \code{"Kurtosis"}.  Here, \code{"Kurtosis"} refers to excess kurtosis (greater than 3), and \code{"PearsonsSkewness"} equals \eqn{\frac{(mean - mode)}{SD}}.  A \code{"Mode"} of 0.5 indicates that the point probabilities at \eqn{x=0} and \eqn{x=1} are tied for highest density; other than this possibility, the LGP is strictly unimodal. 
}
\section{Warning}{
There is a known issue with \code{sLGP()}: when \code{lambda} is negative and \code{theta} is large, the third and fourth moments returned by \code{sLGP()}, with \code{do.numerically=TRUE}, can be quite incorrect due to numerical imprecision.
}
\references{
Consul, P. C.  (1989).  \emph{Generalized Poisson Distributions: Properties and Applications}.  New York:  Marcel Dekker, Inc.

Consul, P. C., & Famoye, F.  (2006).  \emph{Lagrangian Probability Distributions}.  Boston: Birkhauser.

Johnson, N. L., Kemp, A. W., & Kotz, S.  (2005).  \emph{Univariate Discrete Distributions} (3rd. ed.).  Hoboken, NJ: John Wiley & Sons, Inc.
}
\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
LGP.findmax(theta=2, lambda=0.2) #<--No upper support limit
LGP.findmax(theta=2, lambda=-0.2) #<--Upper support limit of 9
LGP.get.nc(theta=2, lambda=0.2)-1==0 #<--TRUE
LGP.get.nc(theta=2, lambda=-0.2)-1 #<--nc differs appreciably from 1
LGP.get.nc(theta=2, lambda=-0.1)-1 #<--nc doesn't differ appreciably from 1
LGPMVP(theta=2, lambda=0.9)
LGPMVP(mu=20, sigma2=2000)
sLGP(theta=2, lambda=0.9)
dLGP(x=0:10,theta=1,lambda=0.1)
dLGP(x=0:10,theta=1,lambda=0)
dLGP(x=0:10,theta=1,lambda=-0.1) #<--Upper support limit of 9
pLGP(q=0:10,theta=1,lambda=0.1)
pLGP(q=0:10,theta=1,lambda=0)
pLGP(q=0:10,theta=1,lambda=-0.1) 
qLGP(p=(0:9)/10,theta=1,lambda=0.1)
qLGP(p=(0:9)/10,theta=1,lambda=0)
qLGP(p=(0:9)/10,theta=1,lambda=-0.1) 
rLGP(n=5,theta=1e12,lambda=-0.0001)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
