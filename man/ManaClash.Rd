\name{ManaClash}
\alias{ManaClash}
\alias{dmanaclash.net}
\alias{dmanaclash.dmg}
\alias{dmanaclash.xyN}
\alias{rmanaclash}

\title{The Mana Clash distributions (just for fun!)}

\description{Density and random-number functions for distributions pertinent to \href{http://gatherer.wizards.com/Pages/Card/Details.aspx?multiverseid=83159}{"Mana Clash"}, a card from the \emph{Magic: The Gathering} trading-card game.  As of 08/29/2014, the official card text read: "You and target opponent each flip a coin. Mana Clash deals 1 damage to each player whose coin comes up tails. Repeat this process until both players' coins come up heads on the same flip."
}
\usage{
dmanaclash.dmg(x,y,N=NULL,pA=0.25,pB=0.25,pC=0.25,pD=0.25,log=FALSE)
dmanaclash.xyN(x,y,N,pA=0.25,pB=0.25,pC=0.25,pD=0.25,log=FALSE)
dmanaclash.net(z,pA=0.25,pB=0.25,pC=0.25,pD=0.25,rel.eps=1e-8,log=FALSE)
rmanaclash(n,pA=0.25,pB=0.25,pC=0.25,pD=0.25,N=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric amount of damage dealt to opponent.}
  \item{y}{Numeric amount of damage dealt to Mana Clash's controller (hereinafter, "you").}
  \item{N}{Numeric.  Number of rounds of coin-tossing that \emph{precede} the last round, when double heads occurs.  That is, \code{N} is the number of rounds of coin-tossing in which at least one player takes damage.}
  \item{z}{Numeric; \emph{net} damage dealt to opponent; negative values are allowed.}
  \item{n}{Integer number of random vectors to generate.}
  \item{pA}{Numeric; probability that both players take damage in a round of coin-tossing (i.e., double tails).}
  \item{pB}{Numeric; probability that you take damage but opponent does not in a round of coin-tossing.}
  \item{pC}{Numeric; probability that opponent takes damage but you do not in a round of coin-tossing.}
  \item{pD}{Numeric; probability that neither player takes damage in a round of coin-tossing (i.e., double heads).}
  \item{log}{Logical; should the natural log of the probability be returned?  Defaults to \code{FALSE}.}
  \item{rel.eps}{Numeric; when computing the sum of an infinite series, how small should the relative change in the sum get before stopping?}
}

\details{
The probability arguments--\code{pA}, \code{pB}, \code{pC}, and \code{pD}--are named as in a two-way contingency table.  They cannot be negative, although \code{rmanaclash} accepts values of zero for \code{pA}, \code{pB}, and \code{pC}.  If they do not sum to 1, they are automatically normalized.  They default to the scenario of two independent fair coins.

Vectors of numeric arguments other than \code{rel.eps} are cycled, whereas only the first element of logical and integer arguments is used.

Function \code{dmanaclash.dmg()} is the bivariate PMF of the amount of damage dealt to opponent and you.  If \code{N = NULL} (default), the probabilities are marginal with respect to the number of rounds of damage-dealing.  Otherwise, the probabilities are conditioned upon the given value of \code{N}.

Function \code{dmanaclash.xyN()} is the trivariate joint PMF of the amount of damage dealt to opponent, the amount dealt to you, and the number of rounds of damage-dealing.

Function \code{dmanaclash.net} is the univariate PMF of the \emph{net} amount of damage dealt to opponent, i.e. damage dealt to opponent minus damage dealt to you.  This distribution has support on the set of integers--including negative values.

Function \code{rmanaclash()} generates random draws from the trivariate joint distribution of \eqn{x}, \eqn{y}, and \eqn{N}; if a non-\code{NULL} value for \code{N} is supplied, the random draws are generated conditionally on that number of damage-dealing rounds.
}
\value{
\code{dmanaclash.dmg()}, \code{dmanaclash.xyN()}, and \code{dmanaclash.net()} all return numeric vectors of probabilities.  \code{rmanaclash()} returns a numeric matrix, with \code{n} rows, and three columns, named "x", "y", and "N".  Each row is a random draw.
}
\references{
\emph{Magic: The Gathering} is a trademark of Wizards of the Coast, LLC, a subsidiary of Hasbro, Inc.
}
\author{Robert M. Kirkpatrick \email{rkirkpatrick2@vcu.edu}}

\section{Derivation}{
Note:  This section is only displayed in the PDF of the package documentation.
\if{latex}{
Let random variables \eqn{X}, \eqn{Y}, and \eqn{N} respectively denote the amount of damage dealt to opponent, the amount of damage dealt to you, and the number of rounds of coin-tossing in which any damage is dealt.  In a given round of coin-tossing, let \eqn{a} denote the probability of damage to both players (i.e., two tails), \eqn{b}, the probability of damage to you but not opponent, \eqn{c}, the probability of damage to opponent but not you, and \eqn{d} the probability of no damage (i.e., two heads, and no more coin-tossing).  Define \eqn{r = a+b+c = 1-d}.

It is obvious that the marginal distribution of \eqn{N} is geometric, with PMF
\deqn{p(n) = dr^n}{p(n) = d*r^n}
Further, given that \eqn{N=n}, the conditional distribution of \eqn{X} is binomial, with PMF
\deqn{p(x | n) = \left(\frac{n!}{(n-x)!x!}\right)\left(\frac{a+c}{r}\right)^x\left(\frac{b}{r}\right)^{n-x}}%{p(x | n) = choose(n,x) * ((a+c)/r)^x * (b/r)^(n-x)}

Now, consider \eqn{Y}, given \eqn{N=n} and \eqn{X=x}.  Of course, \eqn{Y} cannot exceed \eqn{n}.  On any of the \eqn{n} rounds in which the opponent took no damage, you must have taken damage, and you may or may not have taken damage on those rounds in which opponent did.  Thus, \eqn{Y} cannot be less than \eqn{n-x}, and will equal \eqn{n-x} only if exactly one player took damage on all \eqn{n} rounds (i.e., double tails never occurred).  On the \eqn{x} undetermined rounds, you take damage with probability \eqn{\frac{a}{a+c}}{a/(a+c)}, and do not take damage with probability \eqn{\frac{c}{a+c}}{c/(a+c)}.  Thus, the random variable \eqn{Y - (n-x)} follows a binomial distribution, with PMF
\deqn{p(y | x,n) = \left(\frac{x!}{(n-y)!(y-n+x)!}\right)\left(\frac{a}{a+c}\right)^{x+y-n}\left(\frac{c}{a+c}\right)^{n-y}}%{p(y | x,n) = choose(x,y-n+x) * (a/(a+c))^(y-n+x)) * (c/(a+c))^(n-y)}

The joint distribution of \eqn{X} and \eqn{Y}, conditional on \eqn{N=n}, has PMF
\deqn{p(x,y|n) = p(x|n)p(y|x,n)}
%\deqn{= \left(\frac{n!}{(n-x)!x!}\right)\left(\frac{a+c}{r}\right)^x\left(\frac{b}{r}\right)^{n-x} \left(\frac{x!}{(n-y)!(y-n+x)!}\right)\left(\frac{a}{a+c}\right)^{x+y-n}\left(\frac{c}{a+c}\right)^{n-y}}
which is the product of two binomial PMFs.  The joint distribution of \eqn{X}, \eqn{Y}, and \eqn{N} has PMF
\deqn{p(x,y,n) = p(n)p(x|n)p(y|x,n)}
which after some algebraic manipulation can be written as
\deqn{p(x,y,n) = \frac{n!}{(n-x)!(n-y)!(x+y-n)!}b^{n-x}c^{n-y}a^{x+y-n}d}
The joint distribution of \eqn{X} and \eqn{Y}, marginal with respect to \eqn{N}, then has PMF
\deqn{p(x,y) = \frac{da^{x+y}}{b^xc^y}\cdot\sum_{n=\max(x,y)}^{x+y} \left[ \left(\frac{bc}{a}\right)^n \left(\frac{n!}{(n-x)!(n-y)!(x+y-n)!}\right) \right] \ }

Now, define random variable \eqn{Z = X-Y}, i.e. the \emph{net} amount of damage dealt to opponent.  The PMF of \eqn{Z} may be obtained from \eqn{p(x,y)}, the joint PMF of \eqn{X} and \eqn{Y}:
\deqn{p(z) = \sum_{\{(x,y):x-y=z\}}p(x,y) \ }
That is, the probability that \eqn{Z=z} is the sum of the probabilities of all pairs \eqn{(x,y)} such that \eqn{x-y=z}.  For example, the probability that \eqn{Z=-2} is the sum \eqn{p(0,2) + p(1,3) + p(2,4) + ...}
}
}

\examples{
## Same outcome, with and without conditioning on N:
dmanaclash.dmg(x=1,y=1,N=1)
dmanaclash.dmg(x=1,y=1)

## Same damage amounts, with N fixed versus random:
dmanaclash.dmg(x=1,y=1,N=2)
dmanaclash.xyN(x=1,y=1,N=2)

## Net damage distribution is symmetric with defaults:
dmanaclash.net(z=c(-3,-2,-1,0,1,2,3))
## But if coins are biased against opponent...:
dmanaclash.net(z=c(-3,-2,-1,0,1,2,3),pA=0.1,pB=0.1,pC=0.7,pD=0.1)

## Random draws:
rmanaclash(n=10)
rmanaclash(n=10,pA=0.1,pB=0.1,pC=0.7,pD=0.1)
rmanaclash(n=10,N=5)
}
