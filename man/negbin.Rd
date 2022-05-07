\name{negbin}
\alias{dnegbin}
\alias{negbinMVP}

\title{
Helper functions for the (unvariate) negative binomial distribution.
}
\description{
Helper functions for the (unvariate) negative binomial distribution: change-of-parameter, wrapper function for density.
}
\usage{
dnegbin(x,nu,p,mu,log=FALSE)
negbinMVP(nu, p, mu, sigma2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric vector of quantiles.}
  \item{nu}{Numeric; equivalent to argument \code{size} in \code{\link[stats:NegBinomial]{dnbinom}()}, etc.}
  \item{p}{Numeric; equivalent to argument \code{prob} in \code{\link[stats:NegBinomial]{dnbinom}()}, etc.}
  \item{mu}{Numeric vector of mean parameters.}
  \item{sigma2}{"Sigma squared"--numeric vector of variance parameters.}
  \item{log}{Logical; should the natural log of the probability be returned?  Defaults to \code{FALSE}.}
%%     ~~Describe \code{x} here~~
}
\details{
Function \code{dnegbin()} is a wrapper for \code{\link[stats:NegBinomial]{dnbinom}()}.  Two of the three arguments \code{nu}, \code{p}, and \code{mu} must be provided.  Unlike \code{\link[stats:NegBinomial]{dnbinom}()}, \code{dnegbin()} will accept \code{mu} and \code{p}(\code{prob}) with \code{nu}(\code{size}) missing.  In that case, it calculates \code{nu} as \code{mu * p/(1-p)}, and passes \code{nu} and \code{p} to \code{\link[stats:NegBinomial]{dnbinom}()}.

Function \code{negbinMVP()} accepts exactly two of its four arguments, and returns the corresponding values of the other two arguments.  For example, if given values for \code{nu} and \code{p}, it will return the corresponding means (\code{mu}) and variances (\code{sigma2}) of a negative binomial distribution with the given values of \code{nu} and \code{p}.  \code{negbinMVP()} does not strictly enforce the parameter space, but will throw a warning for bad input values.
}

\value{
\code{dnegbin()} returns a numeric vector of probabilities.  \code{negbinMVP()} returns a numeric matrix with two columns, named for the missing arguments in the function call.
}

\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}
}

\seealso{
\code{\link[stats:NegBinomial]{dnbinom}()}
}
\examples{
## These two lines do the same thing:
dnegbin(x=1,nu=2,p=0.5)
dnbinom(x=1,size=2,prob=0.5)

## What is the mean of this distribution?
negbinMVP(nu=2,p=0.5) #<--mu=2

## These two lines also do the same thing:
dnegbin(x=1,nu=2,mu=2)
dnbinom(x=1,size=2,mu=2)

## Parametrize with mu & p:
dnegbin(x=1,mu=2,p=0.5)
## Not run (will throw an error):
## dnbinom(x=1,prob=0.5,mu=2) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
