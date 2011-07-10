\name{runMarkovChainRegenS}
\alias{runMarkovChainRegenS}
\title{
Run a Markov chain with Regeneration on Single Processor
}
\description{
This function provides the capability to run a Markov chain and then split it 
up into independent tours. The user has to provide the minorizing functions 
in order to do so.
}
\usage{
runMarkovChainRegenS(genNextState, genNextStateArgList = NULL, 
initialPoint, regenProbs, regenProbsArgList=NULL, transDens, 
transDensArgList = NULL, smallMeasure, smallMeasureArgList = NULL, 
smallFn, smallFnArgList = NULL, numericParams, nsim, nTours)
}
\arguments{
  \item{genNextState}{The user-specified one-step generation function of the 
  Markov chain. It can have any number of arguments, but one of them must be 
  "x", a vector argument that specifies the current position of the 
  Markov chain.}
  \item{genNextStateArgList}{A named list object, containing the values for the 
  arguments to the user's genNextState() function, other than "x".}
  \item{initialPoint}{The point to start the Markov chain with.}
  \item{regenProbs}{The name of a user-specified function that computes the 
  regeneration probability, given the current and immediate future state. 
  Two of the arguments must be x and y - these denote the current and 
  new states of the chain. If this is specified, then transDens, 
  smallFn and smallMeasure will not be used.}
  \item{regenProbsArgList}{A named list of arguments to the user-specified 
  regenProbs() function. It should contain all arguments except x and y.}
  \item{transDens}{The name of a user-specified function that computes the 
  Markov transition density of the Markov chain. Two of the arguments must be 
  x and y - these vectors denote the current and new states of the chain.}
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
  \item{numericParams}{Not all the variables that are generated are numerical. 
  Some might be categorical. This argument allows one to keep only those of  
  interest.}
  \item{nsim}{Number of iterations of Markov chain to run. As this is one of 2 
  ways of stopping the chian, only one of nsim or nTours must be supplied.}
  \item{nTours}{Number of tours to generate. Only this or nsim must be 
  supplied, not both.}
}
\details{
A Markov chain with transition density \eqn{p(x,y)} can be minorized with 
the following inequality:
\deqn{p(x,y) \ge s(x) \nu(y)}
where \eqn{s} is the small function and \eqn{\nu} is the small measure. With 
these 3 functions, a particular output of the Markov chain can be split into 
independent and identically distributed tours. If the chain is also 
geometrically ergodic, it ensures that we have second moments on the tour 
lengths and the sums of the variables generated within each tour. Thus we can 
use the classical CLT to make inference on the parameters of interest.
}
\value{
If regenerations took place, the function returns a matrix containing the 
information from each tour. Each row in the output matrix corresponds to one 
tour. It will contain the sums of each random variable over that tour, and the 
final column will indicate the length of that tour.

If the minorization is poor, and no regenerations took place, then the output 
matrix consists of ALL the random variables generated. The function 
regenProbsFromRawOutput() can be used on \emph{this} matrix to inspect the regeneration 
probabilities at each iteration.

Remember that if regenProbs argument is provided, then transDens(), 
smallFn() and smallMeasure are not needed.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
 \code{\link{runMarkovChainNoRegenS}}, \code{\link{runMarkovChainRegenS}},
 \code{\link{runMarkovChainRegenP}}, \code{\link{regenProbsFromRawOutput}}
}
\keyword{ htest }
