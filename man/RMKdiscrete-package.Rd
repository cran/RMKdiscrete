\name{RMKdiscrete-package}
\alias{RMKdiscrete-package}
\alias{RMKdiscrete}
\docType{package}

\title{Sundry discrete probability distributions.}

\description{
\emph{RMKdiscrete} implements several univariate and bivariate discrete probability distributions:
  \itemize{
    \item The univariate Lagrangian Poisson distribution has five functions: \code{\link{dLGP}()}, \code{\link{pLGP}()}, \code{\link{qLGP}()}, \code{\link{rLGP}()}, \code{\link{sLGP}()}, \code{\link{LGP.findmax}()}, \code{\link{LGP.get.nc}()}, and \code{\link{LGPMVP}()}.
    \item The bivariate Lagrangian Poisson distribution has three functions: \code{\link{dbiLGP}()}, \code{\link{rbiLGP}()}, and \code{\link{biLGP.logMV}()}.
    \item The bivariate negative binomial distribution has three functions: \code{\link{dbinegbin}()}, \code{\link{rbinegbin}()}, and \code{\link{binegbin.logMV}()}.
    \item Although the \link[stats:NegBinomial]{univariate negative binomial} distribution is implemented in base \R, \emph{RMKdiscrete} provides two helper functions for it: \code{\link{dnegbin}()} and \code{\link{negbinMVP}()}.
  }
Finally, the \code{\link{ManaClash}} distributions are provided just for fun.

The package is presently in an unoptimized but functional "beta" state.  Additional helper functions and distributions, including multivariate distributions in more than two dimensions, are planned for subsequent versions of the package.  Contact the maintainer, \email{rkirkpatrick2@vcu.edu}, with suggestions, bug reports, and feature requests.
}

\details{
\tabular{ll}{
Package: \tab RMKdiscrete\cr
Version: \tab 0.1\cr
Date: \tab 2014/10/17\cr
Depends: \tab R (>= 2.15.0), stats\cr
License: \tab GPL (>= 2)\cr
}}

\author{
Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}
Maintainer: Robert M. Kirkpatrick
}

\keyword{ package }