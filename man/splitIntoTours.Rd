\name{splitIntoTours}
\alias{splitIntoTours}
\title{
Split a Markov chain run into independent tours.
}
\description{
Given output from a Markov chain along with the bell variables, this 
function will split up the chain into independent tours.
}
\usage{
splitIntoTours(outMatrix, bellVars, numericParams)
}
\arguments{
  \item{outMatrix}{Output from a sequential run of a Markov chain, with each 
  row consisting of one iteration of the chain. This is the format that 
  is returned by the function \code{runMarkovChainNoRegenS}, for example.}
  \item{bellVars}{The sampled bell variables that will be used to split the 
  chain. A 1 indicates that the \emph{next} iteration was drawn from the 
  small measure, and hence is independent of the past.}
  \item{numericParams}{Not all the variables that are generated are numerical. 
  Some might be categorical. This argument allows one to keep only those of  
  interest.}
}
\value{
A matrix of tours. It will have one more column than the input argument 
outMatrix. Each row represents one complete tour. The last column denotes the 
length of the tour. The rest of the columns denote the sum of the random 
variables generated during that tour.

If no regenerations were present, NA is returned.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.

  Hobert, J. and Jones, G. and Presnell, B. and Rosenthal, J. (2002) 
  On the applicability of regenerative simulation in Markov chain Monte 
  Carlo. _Biometrika_, *89*, 731.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{internal}
