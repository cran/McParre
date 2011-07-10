\name{runMarkovChainRegenP}
\alias{runMarkovChainRegenP}
\title{
Run a Markov chain with Regeneration on Multiple Processors in Parallel
}
\description{
This function provides the capability to run a Markov chain and then split it 
up into independent tours. The user has to provide the minorizing functions 
in order to do so. This function runs on a cluster.
}
\usage{
runMarkovChainRegenP(genNextState, genNextStateArgList = NULL,
initialPoint, regenProbs, regenProbsArgList=NULL, transDens, 
transDensArgList = NULL, smallMeasure, smallMeasureArgList = NULL, 
smallFn, smallFnArgList = NULL, numericParams, nsim, nTours, prefix, 
pathToLog, outfile="output.rdt")
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
  \item{prefix}{A character string, that prefixes the output log files.}
  \item{pathToLog}{A character string that specifies where to store the 
  output files.}
  \item{outfile}{A character string that specifies the name of the file to 
  which the output matrix should be written.}
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

This function is to be run on a cluster, and uses the MPI paradigm to 
communicate between processors. Each processor will output a log file that 
updates what task it has completed and what it has not. This is useful for 
debugging. This function should be run through R in BATCH mode, so that 
the same function is run on all the nodes. Please see the examples for 
detailed usage. 

Remember that if regenProbs argument is provided, then transDens(), 
smallFn() and smallMeasure are not needed.
}
\value{
No value is returned, but the matrix of tours (or raw output) is saved to the 
specified output file. See \code{\link{runMarkovChainRegenS}} for details on
output matrix.
}
\references{
  Mykland, P. and Tierney, L. and Yu, B. (1995) Regeneration in Markov chain 
  samplers. _Journal of the American Statistical Association_, *90*, 233-241.
}
\author{
Vik Gopal \email{viknesh@stat.ufl.edu}
}
\seealso{
 \code{\link{runMarkovChainRegenS}},
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

### The function is currently defined as
#function(genNextState, genNextStateArgList=NULL, 
#                                 initialPoint, 
#                                 transDens, transDensArgList=NULL,
#                                 smallMeasure, smallMeasureArgList=NULL,
#                                 smallFn, smallFnArgList=NULL, numericParams,
#                                 nsim, nTours) {
#  # Check that exactly one of nsim or nTours is supplied
#  nsimMissing <- missing(nsim)
#  nToursMissing <- missing(nTours)
#  if(nsimMissing & nToursMissing)
#    stop("Either nsim or nTours has to be supplied.\n") else
#  if(!nsimMissing & !nToursMissing)
#    stop("Either nsim or nTours has to be supplied, but not both.\n")
#
#  # check that the number of numeric parameters is known
#  if(missing(numericParams))
#    numericParams <- 1:length(initialPoint)
#
#  # Check the supplied functions and arguments.
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
#
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
## END OF CHECKS ON transDens AND ITS ARGUMENTS. ###############
#  cat("transDens() OK.\n")
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
#  cat("smallMeasure() OK.\n")
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
#    cat("NA's generated by small function.\n")
#
#  # check length of output argument
#  if(length(sX)!=1)
#    stop("Small function should return a real scalar.\n")
## END OF CHECKS ON smallMeasure AND ITS ARGUMENTS. ###############
#  cat("smallMeasure() OK.\n")
#  cat("All checks done.\n\n")
#
############## genOneTour function definition ##################
#genOneTour <- function(startPt, numericParams) {
#  steps <- 1
#  Si <- startPt
#  complete <- 0
#
#   repeat{
#      potentialNewY <- genNextState(x=startPt)
#      regenProb <- smallFn(x=startPt) *
#                   smallMeasure(y=potentialNewY) / 
#                   transDens(x=startPt, y=potentialNewY)
#      complete <- ifelse(runif(1) < regenProb, 1, 0)
#      if (complete==1)
#        break
#      else {
#        Si <- Si + potentialNewY
#        startPt <- potentialNewY
#        steps <- steps + 1
#      }
#    }
#
#  return(list(out=c(Si[numericParams], steps), newStart=potentialNewY))
#  }
################################################################
## if nTours is missing, run the chain for nsim iterations and split that up.
#  if(nToursMissing) {
#    # Run the Markov chain for nsim steps
#    cat("Running Markov chain for", nsim, "iterations .. ", sep=" ")
#    outMatrix <- matrix(0, nrow=nsim, ncol=length(initialPoint))
#    regenProbs <- rep(0, times=nsim)
#    bellVars <- rep(0, times=nsim)
#
#    outMatrix[1,] <- firstPointGenerated
#    for(i in 2:nsim) {
#      outMatrix[i,] <- genNextState(x=outMatrix[i-1,])
#      regenProbs[i-1] <- smallFn(x=outMatrix[i-1,]) * 
#                         smallMeasure(y=outMatrix[i,]) /
#                         transDens(x=outMatrix[i-1,], y=outMatrix[i,])
#    }
#
#    # generate regeneration points
#    bellVars[-nsim] <- ifelse(runif(nsim-1) < regenProbs[-nsim], 1, 0)
#    cat("Done.\n")
#
#    if(sum(bellVars) >= 2) {
#      tourMatrix <- splitIntoTours(outMatrix, bellVars, numericParams)
#      return(tourMatrix)
#    }
#    else {
#      cat("No complete tours found. Returning everything.\n")
#      return(list(outMatrix, regenProbs))
#    }
#  }
#  else { # Run for a specified number of tours
#    # Look for first regeneration point
#    cat("Looking for first regeneration point .. ")
#    findFirstPt <- 0
#    while(findFirstPt == 0) {
#      potentialNewY <- genNextState(x=initialPoint)
#      regenProb <- smallFn(x=initialPoint) *
#                   smallMeasure(y=potentialNewY) /
#                   transDens(x=initialPoint, y=potentialNewY)
#      complete <- ifelse(runif(1) < regenProb, 1, 0)
#      if(complete == 1) {
#        findFirstPt <- 1
#        newStart <- potentialNewY
#      }
#      else {
#        initialPoint <- potentialNewY
#      }
#    }
#    cat("Done.\n")
#
#    # Run the Markov chain for nTours tours
#    cat("Generating", nTours, "tours .. ", sep=" ")
#    tourMatrix <- matrix(0, nrow=nTours, ncol=length(numericParams)+1)
#    for(i in 1:nTours) {
#      tmp <- genOneTour(newStart, numericParams)
#      tourMatrix[i,] <- tmp$out
#      newStart <- tmp$newStart
#    }
#    cat("Done.\n")
#    colnames(tourMatrix) <- c(paste("S", numericParams, sep=""), "N")
#
#    return(tourMatrix)
#
#  }
#
#  }
}
\keyword{ htest }
