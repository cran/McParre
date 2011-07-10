\name{genNextStateOnewayPlain}
\alias{genNextStateOnewayPlain}
\title{
One-step generation function for Oneway Hierarchical Linear Model.
}
\description{
For the oneway hierarchical linear model defined below in the details 
section, this function generates one iteration of the Gibbs sampler to sample 
from the posterior distribution.
}
\usage{
genNextStateOnewayPlain(x, K, M, a1, b1, a2, b2, Y, lam0, mu0)
}
\arguments{
  \item{x}{A vector of length K+3, which is divided as follows. The first 3 
  elements correspond to \eqn{\lambda_\theta}, \eqn{\lambda_e} and 
  \eqn{\mu} respectively. The final K elements correspond to the treatment 
  means.}
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
\details{
This function allows one to run a Gibbs sampler on the following oneway 
hierarchical linear model: (1) \eqn{\lambda_\theta} is drawn from 
Gamma(a1, b1), (2) \eqn{\lambda_e} is drawn from Gamma(a2, b2), (3) \eqn{\mu} 
is drawn from Normal(mu, 1/lam0), (4) \eqn{theta_i} are drawn from 
Normal(\eqn{mu}, \eqn{1/ \lambda_\theta}) and (5) \eqn{y_{ij}} are drawn from 
Normal(\eqn{\theta_i}, \eqn{1/ \lambda_e}). 

The order of the transition that takes place is to generate 
\eqn{\lambda_\theta}, then \eqn{\lambda_e}, then \eqn{\mu} and finally 
the \eqn{\theta_i}'s. For details and for the proof on it's geometric 
ergodicity, please refer to the provided reference below.
}
\value{
A vector of length K+3, which represents the next state of the Gibbs sampler, 
when started at state x.
}
\references{
  Hobert, J. and Geyer, C. (1998) Geometric ergodicity of Gibbs and block 
  Gibbs samplers for a hierarchical random effects model. _Journal of 
  Multivariate Analysis_, *67*, 414-430.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

### The function is currently defined as
#function(x, K, M, a1, b1, a2, b2, Y, lam0, mu0) {
#  y <- rep(0, times=K+3)
#  
#  tmp <- x[-(1:3)]
#  y[1] <- rgamma(1, shape = ( K/2 + a1 ), 
#                 rate = b1 + 0.5 * sum( (tmp - x[3])^2 ) )
#  
#  y[2] <- rgamma(1, shape = ( M*K/2 + a2 ), 
#                 rate = b2 + 0.5 * sum((Y - tmp)^2))
#  
#  y[3] <- rnorm(1, mean = (lam0 * mu0 + K * y[1] * mean(tmp))
#                /(lam0 + K * y[1]), sd =  sqrt(1/(lam0 + K*y[1])))
#  
#  y[-(1:3)] <- rnorm(K, mean = (y[1] * y[3] + M * y[2] * rowMeans(Y))/
#               (y[1] + M * y[2]), sd = sqrt(1/(y[1] + M * y[2])))
#  y
#  }
#}
\keyword{ htest }
