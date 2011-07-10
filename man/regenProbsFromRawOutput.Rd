\name{regenProbsFromRawOutput}
\alias{regenProbsFromRawOutput}
\title{
Compute Regeneration Probabilities for a Markov Chain Run
}
\description{
This is a utility function that allows the user to compute the regeneration 
probabilities for a particular run of a Markov chain. 
}
\usage{
regenProbsFromRawOutput(outMatrix, regenProbs, regenProbsArgList=NULL, 
transDens, transDensArgList = NULL, smallMeasure, 
smallMeasureArgList = NULL, smallFn, smallFnArgList = NULL)
}
\arguments{
  \item{outMatrix}{The output from the Markov chain. It should contain one 
  row for each iteration, and the number of columns should be equal to the 
  number of variables in each state of the chain.}
  \item{regenProbs}{The name of a user-specified function that computes the 
  regeneration probability, given the current and immediate future state. 
  Two of the arguments must be x and y - these denote the current and 
  new states of the chain. If this is specified, then transDens, 
  smallFn and smallMeasure will not be used.}
  \item{regenProbsArgList}{A named list of arguments to the user-specified 
  regenProbs() function. It should contain all arguments except x and y.}
  \item{transDens}{The name of a user-specified function that computes the 
  Markov transition density of the Markov chain. Two of the arguments must be 
  x and y - these denote the current and new states of the chain.}
  \item{transDensArgList}{A named list of arguments to the user-specified 
  transDens() function. It should contain all arguments except x and y.}
  \item{smallMeasure}{The name of a user-specified function that computes the 
  small measure in the minorization of the Markov chain. One of the arguments 
  should be y - this denotes the new state of the chain.}
  \item{smallMeasureArgList}{A named list of arguments to the user-specified 
  smallMeasure() function. It should contain all arguments except y.}
  \item{smallFn}{A user-specified function that computes the small function in 
  the minorization of the Markov chain. One of the arguments must be x - the 
  current state of the chain.}
  \item{smallFnArgList}{A named list of arguments to smallFn() function. It 
  should contain all arguments to smallFn() other than x.}
}
\details{
The regeneration probability is computed as 
\deqn{s(x) \times \nu(y) / p(x,y)}
where \eqn{s} is the small function, \eqn{\nu} is the small measure and 
\eqn{p} is the transition density. This function will typically be run on 
output from \code{runMarkovChainNoRegenS} in order to inspect the quality 
of the minorization before running it in parallel, or before settling on a 
distinguished set.

Remember that if regenProbs argument is provided, then transDens(), 
smallFn() and smallMeasure are not needed.
}
\value{
A vector of regeneration probabilities.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{ htest }
