\name{runMarkovChainRegenIP}
\alias{runMarkovChainRegenIP}
\title{
Interactively Run a Markov chain with Regeneration on Parallel Processors
}
\description{
This function provides the capability to run a Markov chain and then split it 
up into independent tours. The user has to provide the minorizing functions 
in order to do so. This function runs on a cluster. This is different from 
\code{\link{runMarkovChainRegenP}} in that this one runs in an 
interactive Rmpi session.
}
\usage{
runMarkovChainRegenIP(genNextState, genNextStateArgList = NULL, initialPoint, regenProbs, regenProbsArgList = NULL, transDens, transDensArgList = NULL, smallMeasure, smallMeasureArgList = NULL, smallFn, smallFnArgList = NULL, numericParams, nsim, nTours)
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
  \item{nsim}{Number of iterations of Markov chain to run on each processor. 
  As this is one of 2 
  ways of stopping the chian, only one of nsim or nTours must be supplied.}
  \item{nTours}{Number of tours to generate on each processor. Only this or 
  nsim must be supplied, not both.}
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

This function is to be run interactively on a cluster, and uses the MPI 
paradigm to communicate between processors. Before this is run, the slaves 
should already be spawned. Please see the demo3 for detailed usage. 

Remember that if regenProbs argument is provided, then transDens(), 
smallFn() and smallMeasure are not needed.
}
\value{
A matrix containing the generated tours. This tour matrix can be fed into 
\code{\link{computeCI}} in order to obtain Confidence Intervals of interest.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
\code{\link{runMarkovChainRegenP}}
}
\examples{
###---- Should be DIRECTLY executable !! ----
###-- ==>  Define data, use random,
###--	or do  help(data=index)  for the standard data sets.
#
### The function is currently defined as
#function(genNextState, genNextStateArgList=NULL, 
#                       initialPoint, 
#                       regenProbs, regenProbsArgList=NULL,
#                       transDens, transDensArgList=NULL,
#                       smallMeasure, smallMeasureArgList=NULL,
#                       smallFn, smallFnArgList=NULL, numericParams,
#                       nsim, nTours) {
#
#  # Check that slaves are present.
#  commSize <- mpi.comm.size()
#  if(commSize < 2) 
#    stop("No slaves running. Please spawn slaves before calling 
#          runMarkovChainRegenIP().\n") else
#    cat("Found", commSize - 1, "slaves running.\n", sep=" ")
#
#  # Check that exactly one of nsim or nTours is supplied
#  nsimMissing <- missing(nsim)
#  nToursMissing <- missing(nTours)
#  if(nsimMissing & nToursMissing)
#    stop("Either nsim or nTours has to be supplied.\n") else
#  if(!nsimMissing & !nToursMissing)
#    stop("Either nsim or nTours has to be supplied, but not both.\n")
#  if(nsimMissing)
#    mpi.bcast.Robj2slave(nTours) else
#    mpi.bcast.Robj2slave(nsim)
#
#  # Check that either regenProbs, or all of smallFn, smallMeasure and 
#  # transDens are supplied.
#  regenProbsMissing <- missing(regenProbs)
#  if(regenProbsMissing){
#    transDensMissing <- missing(transDens)
#    smallFnMissing <- missing(smallFn)
#    smallMeasureMissing <- missing(smallMeasure)
#    if(transDensMissing || smallFnMissing || smallMeasureMissing)
#      stop("Either regenProbs(), or ALL of transDens(), smallFn() and 
#            smallMeasure() have to be supplied.\n")
#  }
#
#  # check that the number of numeric parameters is known
#  if(missing(numericParams))
#    numericParams <- 1:length(initialPoint)
#  mpi.bcast.Robj2slave(numericParams)
#
## CHECKS ON genNextState AND ITS ARGUMENTS. ######################
#  # check that the genNextState argument is a function.
#  if(class(genNextState)!="function")
#    stop("genNextState() has to be of 'function' class.\n")
#
#  # check that genNextStateArgList is not missing and is a list.
#  if(!is.null(genNextStateArgList)){
#    if(class(genNextStateArgList)!="list") 
#      stop("genNextStateArgList has to be of 'list' class.\n")
#
#    # After this, we should only have to call genNextState(x=...)
#    formals(genNextState) <- c(genNextStateArgList, x=1)
#  }
#
#  # try running one step and stop if there are errors.
#  cat("Checking to see if one-step generation function works .. ")
#  firstPointGenerated <- genNextState(x=initialPoint)
#  if(sum(is.na(firstPointGenerated)) > 0) 
#    stop("NA's generated in genNextState().\n")
#
#  # check length of output argument
#  if(length(firstPointGenerated)!=length(initialPoint))
#    stop("Input/output states for genNextState() do not match in length.\n")
## END OF CHECKS ON genNextState AND ITS ARGUMENTS. ###############
#  cat("genNextState() OK.\n")
#  mpi.bcast.Robj2slave(genNextState)
#
#  if(regenProbsMissing) {
## CHECKS ON transDens AND ITS ARGUMENTS. ######################
#  # check that the genNextState argument is a function.
#  if(class(transDens)!="function")
#    stop("transDens() has to be of 'function' class.\n")
#
#  # check that transDensArgList is not missing and is a list.
#  if(!is.null(transDensArgList)){
#    if(class(transDensArgList)!="list")
#      stop("transDensArgList has to be of 'list' class.\n")
#
#    # After this, we should only have to call genNextState(x=...)
#    formals(transDens) <- c(transDensArgList, x=1, y=1)
#  }
#
#  # try running one step and stop if there are errors.
#  cat("Checking to see if transition density function works .. ")
#  pxy <- transDens(x=initialPoint, y=firstPointGenerated)
#  if(is.na(pxy))
#    stop("NA's generated by transition density function.\n")
#
#  # check length of output argument
#  if(length(pxy)!=1)
#    stop("Transition density function should return a real scalar.\n")
#
## END OF CHECKS ON transDens AND ITS ARGUMENTS. ###############
#  cat("transDens() OK.\n", file=output.con)
#  mpi.bcast.Robj2slave(transDens)
#
## CHECKS ON smallMeasure AND ITS ARGUMENTS. ######################
#  # check that the smallMeasure argument is a function.
#  if(class(smallMeasure)!="function")
#    stop("smallMeasure() has to be of 'function' class.\n")
#
#  # check that smallMeasureArgList is not missing and is a list.
#  if(!is.null(smallMeasureArgList)){
#    if(class(smallMeasureArgList)!="list")
#      stop("smallMeasureArgList has to be of 'list' class.\n")
#
#    # After this, we should only have to call genNextState(x=...)
#    formals(smallMeasure) <- c(smallMeasureArgList, y=1)
#  }
#
#  # try running one step and stop if there are errors.
#  cat("Checking to see if small measure function works .. ")
#  nuY <- smallMeasure(y=firstPointGenerated)
#  if(is.na(nuY))
#    stop("NA's generated by small measure function.\n")
#
#  # check length of output argument
#  if(length(nuY)!=1)
#    stop("Small measure function should return a real scalar.\n")
## END OF CHECKS ON smallMeasure AND ITS ARGUMENTS. ###############
#  cat("smallMeasure() OK.\n", file=output.con)
#  mpi.bcast.Robj2slave(smallMeasure)
#
## CHECKS ON smallFn AND ITS ARGUMENTS. ######################
#  # check that the smallFn argument is a function.
#  if(class(smallFn)!="function")
#    stop("smallFn() has to be of 'function' class.\n")
#
#  # check that smallFnArgList is not missing and is a list.
#  if(!is.null(smallFnArgList)){
#    if(class(smallFnArgList)!="list")
#      stop("smallFnArgList has to be of 'list' class.\n")
#
#    # After this, we should only have to call genNextState(x=...)
#    formals(smallFn) <- c(smallFnArgList, x=1)
#  }
#
#  # try running one step and stop if there are errors.
#  cat("Checking to see if small function works .. ")
#  sX <- smallFn(x=firstPointGenerated)
#  if(is.na(sX))
#    stop("NA's generated by small function.\n")
#
#  # check length of output argument
#  if(length(sX)!=1)
#    stop("Small function should return a real scalar.\n")
## END OF CHECKS ON smallMeasure AND ITS ARGUMENTS. ###############
#  cat("smallFn () OK.\n")
#  mpi.bcast.Robj2slave(smallFn)
#  }
#  else {
## CHECKS ON regenProbs AND ITS ARGUMENTS. ######################
#  # check that the regenProbs argument is a function.
#  if(class(regenProbs)!="function")
#    stop("regenProbs() has to be of 'function' class.\n")
#
#  # When regenProbsArgList is not NULL, check that it is a list.
#  if(!is.null(regenProbsArgList)){
#    if(class(regenProbsArgList)!="list")
#      stop("regenProbsArgList has to be of 'list' class.\n")
#
#    # After this, we should only have to call genNextState(x=...)
#    formals(regenProbs) <- c(regenProbsArgList, x=1, y=1)
#  }
#
#  # try running one step and stop if there are errors.
#  cat("Checking to see if regenProbs function works .. ")
#  firstRegenProb <- regenProbs(x=initialPoint, y=firstPointGenerated)
#  if(is.na(firstRegenProb))
#    stop("NA's generated by regenProbs function with x=", initialPoint, 
#        "and y=", firstPointGenerated, ".\n", sep=" ")
#
#  # check length of output argument
#  if(length(firstRegenProb)!=1)
#    cat("regenProbs function should return a real scalar.\n", file=output.con)
## END OF CHECKS ON regenProbs AND ITS ARGUMENTS. ###############
#  cat("regenProbs() OK.\n")
#  mpi.bcast.Robj2slave(regenProbs)
#  }
#  cat("All checks done.\n\n")
#
#  if(regenProbsMissing) {
#    if(nsimMissing) {
#    A <- mpi.remote.exec(runMarkovChainRegenSslave, genNextState=genNextState, 
#                         initialPoint = initialPoint, 
#                         transDens = transDens, 
#                         smallMeasure = smallMeasure,
#                         smallFn = smallFn, 
#                         numericParams=numericParams, nTours=nTours)
#    } 
#    else {
#    A <- mpi.remote.exec(runMarkovChainRegenSslave, genNextState=genNextState, 
#                         initialPoint = initialPoint, 
#                         transDens = transDens, 
#                         smallMeasure = smallMeasure,
#                         smallFn = smallFn, 
#                         numericParams=numericParams, nsim=nsim)
#
#    }
#  }
#  else {
#    if(nsimMissing) {
#    A <- mpi.remote.exec(runMarkovChainRegenSslave, genNextState=genNextState, 
#                         initialPoint = initialPoint, 
#                         regenProbs = regenProbs,
#                         numericParams=numericParams, nTours=nTours)
#    } 
#    else {
#    A <- mpi.remote.exec(runMarkovChainRegenSslave, genNextState=genNextState, 
#                         initialPoint = initialPoint, 
#                         regenProbs = regenProbs,
#                         numericParams=numericParams, nsim=nsim)
#    }
#  }
#
#  if(nsimMissing) {
#    tourMatrix <- NULL
#    for(ii in 1:(commSize-1)) {
#      tmpTours <- nrow(A[[ii]])
#      cat(tmpTours, "tours from slave", ii, ".\n", sep=" ")
#      tourMatrix <- rbind(tourMatrix, A[[ii]])
#    }
#    return(tourMatrix)
#  }
#  else {
#    tourMatrix <- NULL
#    for(ii in 1:(commSize - 1)) {
#      if(!("N" \%in\% colnames(A[[ii]]))) {
#        cat("No tours from slave", ii, ".\n", sep=" ")
#        next
#      }
#      tmpTours <- nrow(A[[ii]])
#      cat(tmpTours, "tours from slave", ii, ".\n", sep=" ")
#      tourMatrix <- rbind(tourMatrix, A[[ii]])
#    }
#    if(is.null(tourMatrix)) {
#      cat("Returning raw output from slave 1.\n")
#      return(A[[1]])
#    }
#    else 
#      return(tourMatrix)
#  }
#
#  }
}
\keyword{ htest }
