\name{smallMeasureHLMM}
\alias{smallMeasureHLMM}
\title{
Compute Small Measure in Minorization of Hierarchical Linear Mixed Model.
}
\description{
For the hierarchical linear mixed model minorized with a distinguished point 
and set as described in the reference below, this function computes the 
value of the small measure, which is needed when computing regeneration
probabilities.
}
\usage{
smallMeasureHLMM(xTil, y, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n, lamRInterval, lamDInterval)
}
\arguments{
  \item{xTil}{The distinguished point.}
  \item{y}{A vector of length p+q+2, which is divided as follows. The first 
  two elements correspond to \eqn{\lambda_D} and \eqn{\lambda_R} 
  respectively. The next q elements represent the u vector, while the last 
  p elements make up the \eqn{\beta} vector. This is the state that will 
  be tested to see if it is in the distinguished set.}
  \item{Y}{Vector of length n, contains the observations.}
  \item{X}{Design matrix for ``fixed'' effects. It must have n rows and 
  p columns. It must have full column rank (check that X^t X is 
  invertible).}
  \item{Z}{Design matrix for random effects. It must have n rows and 
  q columns.}
  \item{r1}{Shape parameter for prior on \eqn{\lambda_R}.}
  \item{r2}{Rate parameter for prior on \eqn{\lambda_R}.}
  \item{d1}{Shape parameter for prior on \eqn{\lambda_D}.}
  \item{d2}{Rate parameter for prior on \eqn{\lambda_D}.}
  \item{beta0}{Prior mean for \eqn{\beta} vector.}
  \item{Bmat}{Inverse of the covariance matrix for the prior distribution 
  on \eqn{\beta}.}
  \item{p}{Number of columns of X-matrix = Length of \eqn{\beta}
  vector.}
  \item{q}{Number of columns of Z-matrix = Length of u-vector.}
  \item{n}{Number of observations.}
  \item{lamRInterval}{The distinguished set for \eqn{\lambda_R}.}
  \item{lamDInterval}{The distinguished set for \eqn{\lambda_D}.}
}
\value{
A real value, denoting the value of the small measure for 
this minorization.  }
\references{
  Hobert, J. and Jones, G. and Robert, C. (2006) Using a Markov chain to 
  construct a tractable aproximation of an intractable probability 
  distribution. _Scandinavian Journal of Statistics_, *33*, 37-51.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{smallFnHLMM}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(xTil, y, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, 
#                             q, n, lamRInterval, lamDInterval) {
#  tmp1 <- findInterval(y[1], lamDInterval)
#  tmp2 <- findInterval(y[2], lamRInterval)
#  if((tmp1!=1)||(tmp2!=1))
#    return(0)
#
#  dgamma(y[1], q/2 + d1, rate = d2 + 0.5 * drop(crossprod(xTil[3:(2+q)]))) * 
#  dgamma(y[2], n/2 + r1, rate = r2 + 0.5 * 
#                 drop(crossprod(Y - X \%*\% xTil[(2+q+1):(q+p+2)] - 
#                 Z \%*\% xTil[3:(2+q)])))
#  }
}
\keyword{ htest }
