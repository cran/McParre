\name{genNextStateHLMM}
\alias{genNextStateHLMM}
\title{
One-step generation function for Hierarchical Linear Mixed Model.
}
\description{
For the hierarchical linear mixed model defined below in the details 
section, this function generates one iteration of the Gibbs sampler to sample 
from the posterior distribution.
}
\usage{
genNextStateHLMM(x, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n)
}
\arguments{
  \item{x}{A vector of length p+q+2, which is divided as follows. The first 
  two elements correspond to \eqn{\lambda_D} and \eqn{\lambda_R} 
  respectively. The next q elements represent the u vector, while the last 
  p elements make up the \eqn{\beta} vector.}
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
}
\details{
This function allows one to run a Gibbs sampler on the
hierarchical linear mixed model described in the reference 
below.
}
\value{
A vector of length p+q+2, which represents the next state of the Gibbs sampler, 
when started at state x.
}
\references{
  Hobert, J. and Jones, G. and Robert, C. (2006) Using a Markov chain to 
  construct a tractable aproximation of an intractable probability 
  distribution. _Scandinavian Journal of Statistics_, *33*, 37-51.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(x, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n) {
#  require(MASS)
#  y <- rep(0, q+p+2)
#
#  # generate lamD variable - recall that x is 
#  # arranged as (lamD, lamR, u, beta)
#  y[1] <- rgamma(1, q/2 + d1, rate=d2 + 0.5 * drop(crossprod(x[3:(2+q)])))
#
#  # generate lamR
#  y[2] <- rgamma(1, n/2 + r1, rate = r2 + 0.5 * 
#                 drop(crossprod(Y - X \%*\% x[(2+q+1):(q+p+2)] - 
#                 Z \%*\% x[3:(2+q)])))
#
#  # generate xi
#  Rmat <- y[2] * diag(n)
#  Dmat <- y[1] * diag(q)
#  tmpSigma1 <- cbind(t(Z) \%*\% Rmat \%*\% Z + Dmat, t(Z) \%*\% Rmat \%*\% X)
#  tmpSigma2 <- cbind(t(X) \%*\% Rmat \%*\% Z, t(X) \%*\% Rmat \%*\% X + Bmat)
#  covMat <- solve(rbind(tmpSigma1, tmpSigma2))
#
#  meanVec <- drop(covMat \%*\% rbind(t(Z) \%*\% Rmat \%*\% Y,  
#                              t(X) \%*\% Rmat \%*\% Y + Bmat \%*\% beta0))
#
#  y[-(1:2)] <- mvrnorm(1, mu=meanVec, Sigma=covMat)
#
#  y
#  }
}
\keyword{ htest }
