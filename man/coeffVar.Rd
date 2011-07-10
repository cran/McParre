\name{coeffVar}
\alias{coeffVar}
\title{
Compute Coefficient of Variation of Generated Tours
}
\description{
The Coefficient of Variation (CV) of the generated tours has been mentioned 
as a possible way to see if sufficient tours have been generated before the 
CLT is applied. This function will compute the estimate of the CV from the 
matrix of tours (the output of \code{runMarkovChainRegenS}, for example).
}
\usage{
coeffVar(outObj)
}
\arguments{
  \item{outObj}{A matrix of tours returned from \code{runMarkovChainRegenP} or 
  \code{runMarkovChainRegenS}. In fact, any matrix can be submitted to the 
  function as long as it has a column named "N". This is the only column that 
  the function utilises.}
}
\value{
  A scalar that is the estimate of the CV of the tours that were generated.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{ htest }
