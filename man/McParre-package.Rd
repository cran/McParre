\name{McParre-package}
\alias{McParre-package}
\alias{McParre}
\docType{package}
\title{
Parallel Regenerative Markov Chains
}
\description{
Allows a user to run a regenerative Markov chain (MC) in parallel. If a MC 
can be shown to be geometricaly ergodic, it can be ``split'' into 
i.i.d tours, to which a Central Limit Theorem can be applied. This package 
allows a user to run those tours in parallel, on multiple processors.
}
\details{
\tabular{ll}{
Package: \tab McParre\cr
Type: \tab Package\cr
Version: \tab 0.20\cr
Date: \tab 2011-07-10\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
The user has to provide 4 basic functions - the one-step generation function, 
the Markov transition density, the small function, and the small measure. The 
main functions in the package is \code{runMarkovChainRegenP}. This will 
run a specified number of tours on each processor, or run independent 
chains on seperate processors and then split them up into independent 
tours. However, the user should first run \code{runMarkovChainNoRegenS} and 
\code{runMarkovChainRegenS} on a single processor to ensure that the functions 
are working correctly.

As examples, we provide functions to run a oneway hierarchical linear model 
and a mixed effects hierarchical linear model.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
Maintainer: Vik Gopal <viknesh@stat.ufl.edu>
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.

  Hobert, J. and Jones, G. and Presnell, B. and Rosenthal, J. (2002) 
  On the applicability of regenerative simulation in Markov chain Monte 
  Carlo. _Biometrika_, *89*, 731.
}
\keyword{ package }
\seealso{
 \code{\link{runMarkovChainNoRegenS}}, \code{\link{runMarkovChainRegenS}},
 \code{\link{runMarkovChainRegenP}}
}
