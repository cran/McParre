\name{regenProbsOnewayDiffuse}
\alias{regenProbsOnewayDiffuse}
\title{
Compute Regeneration Probability Given Current and Future States
}
\description{
For the oneway hierarchical linear model minorized with a distinguished point 
and set as described in the reference below, this function computes the value 
of the small measure, which is needed when computing regeneration
probabilities. This models assumes improper priors.
}
\usage{
regenProbsOnewayDiffuse(x, y, xTil, Dmat, a, b, Y, q, miVec, trtMeans)
}
\arguments{
  \item{x}{ A vector of length q+3, which is divided as follows. The first 2
  elements correspond to \eqn{\sigma^2_\theta} and 
  \eqn{\sigma^2_e} respectively. The final q+1 elements correspond 
  to the sub-population \eqn{\theta_i} and \eqn{\mu}. }
  \item{y}{ The immediate future state. }
  \item{xTil}{ The distinguished point. The length of this vector should be 
  same as x and y.}
  \item{Dmat}{ The distinguished set to check if \eqn{\sigma^2_\theta} and 
  \eqn{\sigma^2_e} are in.}
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
\value{
A real value, denoting the regeneration probability for this minorization when 
moving from state x to y.  }
\references{
  Tan, A. and Hobert, J.P. (2009) Block Gibbs sampling for Bayesian random 
  effects models with improper priors: convergence and regeneration.
  _Journal of Computational and Graphical Statistics_, *18*, 861-878.

  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{genNextStateOnewayDiffuse}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
#function(x, y, xiTil, Dmat, a, b, Y, q, miVec, trtMeans) {
#
#  if(length(xiTil) != length(x) - 3)
#    stop("Length of distinguished point is incorrect.\n")
#
#  # check if y is in distinguished set
#  tmp <- apply(cbind(y[1:2], Dmat), 1, function(x)
#               {findInterval(x[1], x[-1]) == 1})
#  if(sum(tmp) < 2) return(0)
#
#  # remember, order is sigT, sigE, mu, thetaVec
#  w1 <- drop(crossprod(x[-(1:3)] - x[3]))
#  w2 <- drop(crossprod(miVec, (trtMeans - x[-(1:3)])^2))
#
#  w1Star <- drop(crossprod(xiTil[-1] - xiTil[1]))
#  w2Star <- drop(crossprod(miVec, (trtMeans - xiTil[-1])^2))
#
#  if(w1 > w1Star) 
#    sigT <- Dmat[1,1] else 
#    sigT <- Dmat[1,2]
#  if(w2 > w2Star) 
#    sigE <- Dmat[2,1] else
#    sigE <- Dmat[2,2]
#
#  logProb <- 0.5 * ( (w1 - w1Star)*(1/y[1] - 1/sigT) + 
#                     (w2 - w2Star)*(1/y[2] - 1/sigE)  )
#
#  exp(logProb)
#  }
}
\keyword{ htest }
