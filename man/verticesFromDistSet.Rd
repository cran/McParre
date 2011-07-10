\name{verticesFromDistSet}
\alias{verticesFromDistSet}
\alias{smallFnComponent}
\alias{checkIfInDistSet}
\title{
Utility functions for Topic Model chain
}
\description{
These functions are internal to regenProbsTopicModels.
}
\usage{
verticesFromDistSet(distSet)
smallFnComponent(coeff, vertices)
checkIfInDistSet(vec1, distSet)
}
\arguments{
  \item{distSet}{Distinguished set for theta or beta.  }
  \item{coeff}{Coefficients in computing the small function.}
  \item{vertices}{Vertices over which the infimum is taken.}
  \item{vec1}{The vector to be checked to see if it is in the distinguished 
  set.}
}
\value{
 \code{verticesFromDistSet} returns a matrix. 
 \code{checkIfInDistSet} returns a TRUE or FALSE.
 \code{smallFnComponent} returns a real value.
}
\keyword{ internal }
