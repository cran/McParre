\name{genNextStateOnewayDiffuse}
\alias{genNextStateOnewayDiffuse}
\title{
One-step generation function for Oneway Model with Diffuse Priors
}
\description{
For the oneway linear model defined below in the details 
section, this function generates one iteration of the Gibbs sampler to sample 
from the posterior distribution. This implementation allows for unbalanced 
data.
}
\usage{
genNextStateOnewayDiffuse(x, a, b, Y, q, miVec, trtMeans)
}
\arguments{
  \item{x}{ A vector of length q+3, which is divided as follows. The first 2
  elements correspond to \eqn{\sigma^2\theta} and 
  \eqn{\sigma^2} respectively. The final q+1 elements correspond 
  to the sub-population \eqn{\theta_i} and \eqn{\mu}. }
  \item{a}{ The parameter for the uninformative prior on \eqn{\sigma^2_\theta}}
  \item{b}{The parameter for the uninformative prior on \eqn{\sigma^2_e}}
  \item{Y}{The observed data. This should be a matrix of with K rows and M 
  columns. Missing values should be replaced with an NA.}
  \item{q}{ The number of treatments. }
  \item{miVec}{ The number of observation under each treatment. This should 
  be a vector of length q.}
  \item{trtMeans}{ The mean of the observations under each treatment. This 
  should be a vector of length q. Each element will be the corresponding 
  row mean of Y.}
}
\details{
This function allows one to run a Gibbs sampler on the following oneway 
hierarchical linear model: (1) \eqn{\sigma^2_\theta \propto (\sigma^2)^(-a-1)} 
(2) \eqn{\sigma^2_e \propto (\sigma^2)^(-b-1)}, 
(3) \eqn{\mu \propto 1} 
(4) \eqn{theta_i} are drawn from 
Normal(\eqn{mu}, \eqn{\sigma^2_\theta}) and (5) \eqn{y_{ij}} are drawn from 
Normal(\eqn{\theta_i}, \eqn{\sigma^2_e}). 

The order of the transition that takes place is to generate 
\eqn{\sigma^2_\theta}, then \eqn{\sigma^2_e}, and then the \eqn{\mu} and 
the \eqn{\theta_i}'s as a block. For details and for the proof on it's 
geometric ergodicity, please refer to the provided reference below.
}
\value{
A vector of length q+3, which represents the next state of the Gibbs sampler, 
when started at state x.
}
\references{
  Tan, A. and Hobert, J.P. (2009) Block Gibbs sampling for Bayesian random 
  effects models with improper priors: convergence and regeneration.
  _Journal of Computational and Graphical Statistics_, *18*, 861-878.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{regenProbsOnewayDiffuse}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(x, a, b, Y, q, miVec, trtMeans) {
#
#  newPt <- rep(0, length(x))
#
#  # generate sigT
#  newPt[1] <- 1/rgamma(1, q/2 + a, rate=drop(0.5*crossprod(x[-(1:3)] - x[3])) )
#
#  # generate sigE
#  newPt[2] <- 1/rgamma(1, sum(miVec)/2 + b, 
#             rate=sum( (Y - x[-(1:3)])^2, na.rm=TRUE))
#
#  # set up meanVec
#  tmp1 <- miVec/(newPt[2] + miVec*newPt[1])
#  tmp2 <- drop(crossprod(tmp1, trtMeans)/sum(tmp1))
#  meanVec <- c(tmp2, newPt[2]/(newPt[2] + miVec*newPt[1]) * tmp2 + 
#                     tmp1 * newPt[1] * trtMeans)
#
#  tmp2 <- c(1, newPt[2] / ( newPt[2] + miVec*newPt[1] ))
#  sigMat <- 1/sum(tmp1) * as.matrix(tmp2) \%*\% t(as.matrix(tmp2))
#  diag(sigMat) <- c(1/sum(tmp1), tmp2[-1] * (newPt[1] + tmp2[-1]/sum(tmp1)) )
#
#  newPt[-(1:2)] <- mvrnorm(1, mu=meanVec, Sigma=sigMat)
#  newPt 
#  }
}
\keyword{ htest }
