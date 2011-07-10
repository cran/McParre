\name{transDensOnewayPlain}
\alias{transDensOnewayPlain}
\title{
Compute Markov transition density for Oneway Hierarchical Linear Model.
}
\description{
For the oneway hierarchical linear model defined below in the details 
section, this function computes the transition density for the move from state 
x to y.
}
\usage{
transDensOnewayPlain(x, y, K, M, a1, b1, a2, b2, Y, lam0, mu0)
}
\arguments{
  \item{x}{A vector of length K+3, which is divided as follows. The first 3 
  elements correspond to \eqn{\lambda_\theta}, \eqn{\lambda_e} and 
  \eqn{\mu} respectively. The final K elements correspond to the treatment 
  means. This is the starting state.}
  \item{y}{Same as x, but this is the \emph{next} step.}
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
This function computes the transition density for 
a Gibbs sampler on the following oneway 
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
A real value, representing the product of the conditional densities in the 
above Gibbs sampler.
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

## The function is currently defined as
#function(x, y, K, M, a1, b1, a2, b2, Y, lam0, mu0) {
#    tmp <- x[-(1:3)]
#
#    p1 <- dgamma(y[1], shape = ( K/2 + a1 ), 
#                 rate = b1 + 0.5 * sum( (tmp - x[3])^2 )) *
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% __ONLY ONE__ keyword per line
