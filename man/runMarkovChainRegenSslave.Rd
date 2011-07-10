\name{runMarkovChainRegenSslave}
\alias{runMarkovChainRegenSslave}
\title{
Run a Regenerative Markov chain on Slave Nodes.
}
\description{
This function is called by \code{\link{runMarkovChainRegenIP}} and is not meant
for calling by the user.
}
\usage{
runMarkovChainRegenSslave(genNextState, initialPoint, regenProbs, transDens, smallMeasure, smallFn, numericParams, nsim, nTours)
}
\arguments{
  \item{genNextState}{The user-specified one-step generation function of the 
  Markov chain. It can have any number of arguments, but one of them must be 
  "x", a vector argument that specifies the current position of the 
  Markov chain.}
  \item{initialPoint}{The point to start the Markov chain with.}
  \item{regenProbs}{The name of a user-specified function that computes the 
  regeneration probability, given the current and immediate future state. 
  Two of the arguments must be x and y - these denote the current and 
  new states of the chain. If this is specified, then transDens, 
  smallFn and smallMeasure will not be used.}
  \item{transDens}{The name of a user-specified function that computes the 
  Markov transition density of the Markov chain. Two of the arguments must be 
  x and y - these vectors denote the current and new states of the chain.}
  \item{smallMeasure}{The name of a user-specified function that computes the 
  small measure in the minorization of the Markov chain. One of the arguments 
  should be y - this denotes the new state of the chain.}
  \item{smallFn}{A user-specified function that computes the small function in 
  the minorization of the Markov chain. One of the arguments must be x - the 
  current state of the chain.}
  \item{numericParams}{Not all the variables that are generated are numerical. 
  Some might be categorical. This argument allows one to keep only those of  
  interest.}
  \item{nsim}{Number of iterations of Markov chain to run on each processor. 
  As this is one of 2 
  ways of stopping the chian, only one of nsim or nTours must be supplied.}
  \item{nTours}{Number of tours to generate on each processor. Only this or 
  nsim must be supplied, not both.}
}
\value{
A matrix containing the tours from individual processors.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{ internal }
