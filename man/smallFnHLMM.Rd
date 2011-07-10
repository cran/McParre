\name{smallFnHLMM}
\alias{smallFnHLMM}
\title{
Compute Small Function in Minorization of Hierarchical Linear Mixed Model.
}
\description{
For the hierarchical linear mixed model minorized with a distinguished point 
and set as described in the reference below, this function computes the 
value of the small function, which is needed when computing regeneration
probabilities.
}
\usage{
smallFnHLMM(x, xTil, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n, lamRInterval, lamDInterval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector of length p+q+2, which is divided as follows. The first 
  two elements correspond to \eqn{\lambda_D} and \eqn{\lambda_R} 
  respectively. The next q elements represent the u vector, while the last 
  p elements make up the \eqn{\beta} vector. This is the starting state.}
  \item{xTil}{The distinguished point.}
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
A real value, denoting the value of the small function for 
this minorization.}
\references{
  Hobert, J. and Jones, G. and Robert, C. (2006) Using a Markov chain to 
  construct a tractable aproximation of an intractable probability 
  distribution. _Scandinavian Journal of Statistics_, *33*, 37-51.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{smallMeasureHLMM}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(x, xTil, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, 
#                        n, lamRInterval, lamDInterval) {
#  
#  ux <- x[3:(2+q)]
#  uxTil <- xTil[3:(2+q)]
#
#  fac1 <- drop(crossprod(ux) - crossprod(uxTil))
#  if (fac1 > 0)
#    rat1 <- lamDInterval[2] else 
#    rat1 <- lamDInterval[1]
#  tmp1 <- (d2 + 0.5 * drop(crossprod(ux)))/(d2 + 0.5 * drop(crossprod(uxTil)))
#  tmp1 <- tmp1^(q/2 + d1) * exp(-0.5*fac1*rat1)
#
#  vx <- Y - X \%*\% x[(2+q+1):(q+p+2)] - Z \%*\% x[3:(2+q)]
#  vxTil <- Y - X \%*\% xTil[(2+q+1):(q+p+2)] - Z \%*\% xTil[3:(2+q)]
#
#  fac2 <- drop(crossprod(vx) - crossprod(vxTil))
#  if (fac2 > 0)
#    rat2 <- lamRInterval[2] else 
#    rat2 <- lamRInterval[1]
#  tmp2 <- (r2 + 0.5 * drop(crossprod(vx)))/(r2 + 0.5 * drop(crossprod(vxTil)))
#  tmp2 <- tmp2^(n/2 + r1) * exp(-0.5*fac2*rat2)
#
#  tmp1 * tmp2
#  }
}
\keyword{ htest }
