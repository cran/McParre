\name{srq}
\alias{srq}
\title{
Plot Scaled Regeneration Quantile Plot of Generated Tours
}
\description{
The scaled regeneration quantile (SRQ) plot has 
been mentioned as a possible way to 
see if sufficient tours have been generated before the CLT is applied. 
This function will plot the SRQ plot.
}
\usage{
srq(outputObj, subset, ...)
}
\arguments{
  \item{outputObj}{A matrix of tours returned from \code{runMarkovChainRegenP} 
  or 
  \code{runMarkovChainRegenS}. In fact, any matrix can be submitted to the 
  function as long as it has a column named "N". This is the only column that 
  the function utilises.}
  \item{subset}{ A vector of indices specifying which tours to use for the 
  plot. The reason for this is that typically we have thousands of tours and 
  so a plot with all of them will not be clear.  }
  \item{\dots}{ Extra arguments to be passed to the plot() function. }
}
\details{
  The plot comes from \eqn{T_i/T_R} versus \eqn{i/R} where \eqn{T_i} is the 
  number of iterations to reach the i-th tour, and there are \eqn{R} tours 
  in total.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\keyword{ htest }
