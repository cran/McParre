\name{smallMeasureOnewayPlain}
\alias{smallMeasureOnewayPlain}
\title{
Compute Small Measure in Minorization of Oneway Hierarchical Linear Model.
}
\description{
For the oneway hierarchical linear model minorized with a distinguished point 
and set as described in the reference below, this function computes the value 
of the small measure, which is needed when computing regeneration
probabilities.
}
\usage{
smallMeasureOnewayPlain(xTil, y, D, K, M, a1, b1, a2, b2, Y, lam0, mu0)
}
\arguments{
  \item{xTil}{The distinguished point.}
  \item{y}{A vector of length K+3, which is divided as follows. The first 3 
  elements correspond to \eqn{\lambda_\theta}, \eqn{\lambda_e} and 
  \eqn{\mu} respectively. The final K elements correspond to the treatment 
  means. This is the state that will be tested to see if it is in the 
  distinguished set D}
  \item{D}{The distinguished set. This should be a matrix with 3 rows and 2 
  columns. The rows correspond to intervals for \eqn{\lambda_\theta}, 
  \eqn{\lambda_e} and \eqn{\mu} respectively. If y[1:3] fall within these
  intervals, then the value of the small measure is equal to the 
  transition density from xTil to y. Otherwise, it is 0.} 
  \item{K}{The number of treatments.}
  \item{M}{The number of observations per treatment.}
  \item{a1}{Parameter for prior on \eqn{\lambda_\theta}.}
  \item{b1}{Parameter for prior on \eqn{\lambda_\theta}.}
  \item{a2}{Parameter for prior on \eqn{\lambda_e}.}
  \item{b2}{Parameter for prior on \eqn{\lambda_e}.}
  \item{Y}{The observed data. This should be a matrix of with K rows and M 
  columns.}
  \item{lam0}{Parameter for prior on \eqn{\mu}.}
  \item{mu0}{Parameter for prior on \eqn{\mu}.}
}
\value{
A real value, denoting the value of the small measure for this minorization.  }
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.

  Hobert, J. and Geyer, C. (1998) Geometric ergodicity of Gibbs and block 
  Gibbs samplers for a hierarchical random effects model. _Journal of 
  Multivariate Analysis_, *67*, 414-430.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{transDensOnewayPlain}}, \code{\link{genNextStateOnewayPlain}},
\code{\link{smallFnOnewayPlain}}
}
\examples{
###---- Should be DIRECTLY executable !! ----
###-- ==>  Define data, use random,
###--	or do  help(data=index)  for the standard data sets.
#
### The function is currently defined as
#function(y, xTil, D, K, M, a1, b1, 
#                                  a2, b2, Y, lam0, mu0) {
#
#    tmp <- apply(cbind(y[1:3], D), 1, function(x) 
#                 {findInterval(x[1], x[-1]) == 1})
#    if(sum(tmp) < 3) return(0)
#
#    tmp <- xTil[-(1:3)]
#
#    p1 <- dgamma(y[1], shape = ( K/2 + a1 ), 
#                 rate = b1 + 0.5 * sum( (tmp - xTil[3])^2 )) *
#          dgamma(y[2], shape = ( M*K/2 + a2 ), 
#                   rate = b2 + 0.5 * sum( (Y - tmp)^2 ) ) *
#          dnorm(y[3], mean = (lam0*mu0 + K*y[1]*mean(tmp))/
#                  (lam0 + K*y[1]), sd =  sqrt(1/(lam0 + K*y[1])))
#    p2 <- prod(dnorm(y[-(1:3)], mean = 
#                     (y[1]*y[3] + M*y[2]*rowMeans(Y))/(y[1] + M*y[2]), 
#                     sd = sqrt(1/(y[1] + M*y[2]))))
#
#  p1 * p2
#  }
}
\keyword{ htest }% __ONLY ONE__ keyword per line
