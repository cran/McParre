\name{computeCI}
\alias{computeCI}
\title{
Compute Confidence Intervals on Parameters of Interest
}
\description{
Once tours have been generated, they will be stored in a matrix, with the 
information in one row corresponding to one tour. The final column contains 
the length of the tour, and the other columns contain the sum of the 
random variable generated during that tour. This function computes 
the Confidence Interval at the specified alpha-level for parameters of 
interest.
}
\usage{
computeCI(outObj, alpha = 0.05, paramsOfInterest)
}
\arguments{
  \item{outObj}{This should be a matrix of tours, as returned by 
  \code{runMarkovChainRegenP} or \code{runMarkovChainRegenS}.}
  \item{alpha}{The alpha level of the Confidence Interval to be computed.}
  \item{paramsOfInterest}{Not all random variables in a Markov chain are 
  numeric, even if they have to be coded as such here. For example, some might 
  be categorical variables. As it would not make sense for Confidence 
  Intervals to be computed on those, this argument allows the user to 
  specify which random variables in the matrix to make inference on. If this 
  argument is left unspecified, then the function will compute Confidence 
  Intervals for all variables in the tour matrix.}
}
\value{
The function returns a matrix of 3 columns, and with the number of rows equal 
to the length of paramsOfInterest. Each row corresponds to a parameter for 
which a Confidence Interval has been computed. The first and third columns are 
the lower and upper end-points of the Confidence Interval, while the middle 
column is the point estimate of the parameter.
}
\references{
  Hobert, J. and Jones, G. and Presnell, B. and Rosenthal, J. (2002) 
  On the applicability of regenerative simulation in Markov chain Monte 
  Carlo. _Biometrika_, *89*, 731.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{ htest }
