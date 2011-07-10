\name{runMarkovChainNoRegenS}
\alias{runMarkovChainNoRegenS}
\title{
Run a Markov chain Without Regeneration on Single Processor
}
\description{
This function allows a user to run a vanilla Markov chain for a specified 
number of iterations. The input arguments should be the one-step generation 
function of the desired Markov chain.
}
\usage{
runMarkovChainNoRegenS(genNextState, genNextStateArgList = NULL, initialPoint, nsim)
}
\arguments{
  \item{genNextState}{The one-step generation function of the Markov chain. It 
  can have any number of arguments, but one of them should be "x", which 
  specifies the current position of the Markov chain.}
  \item{genNextStateArgList}{A named list object, containing the values for the 
  arguments to the user's genNextState() function, other than "x".}
  \item{initialPoint}{The start point of the Markov chain.}
  \item{nsim}{The number of iterations of the Markov chain to run.}
}
\details{
Typically, the user's function for one step generation would have many 
arguments, e.g. hyperprior settings, data, etc. These would not change during 
the course of one run of a Markov chain. The second argument above allows 
these settings to be passed to genNextState(). Within runMarkovChainNoRegenS(),
the formal arguments to genNextState() are fixed, this making it only a 
function of x for the duration of the run.
}
\value{A matrix consisting of nsim rows and number of columns equal to the 
length of the initial point. Each row represents one iteration of the chain.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\note{
This a serial function - meant to be run on a single processor. The output is 
meant to be used for determining the distinguished point/set for a particular 
chain.
}
\keyword{ htest }
