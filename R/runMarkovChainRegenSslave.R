runMarkovChainRegenSslave <-
function(genNextState, initialPoint, regenProbs, transDens, smallMeasure, 
                       smallFn, numericParams, nsim, nTours) {
  nsimMissing <- missing(nsim)
  nToursMissing <- missing(nTours)
  regenProbsMissing <- missing(regenProbs)

############# genOneTour function definition ##################
genOneTour <- function(startPt, numericParams) {
  steps <- 1
  Si <- startPt
  complete <- 0

   repeat{
      potentialNewY <- genNextState(x=startPt)
      if(regenProbsMissing) {
        regenProb <- smallFn(x=startPt) *
                     smallMeasure(y=potentialNewY) /
                     transDens(x=startPt, y=potentialNewY)
      }
      else {
        regenProb <- regenProbs(x=startPt, y=potentialNewY)
      }
      complete <- ifelse(runif(1) < regenProb, 1, 0)
      if (complete==1)
        break
      else {
        Si <- Si + potentialNewY
        startPt <- potentialNewY
        steps <- steps + 1
      }
    }

  return(list(out=c(Si[numericParams], steps), newStart=potentialNewY))
}
###############################################################
# if nTours is missing, run the chain for nsim iterations and split that up.
  if(nToursMissing) {
    # Run the Markov chain for nsim steps
    # cat("Running Markov chain for", nsim, "iterations .. ", sep=" ")
    outMatrix <- matrix(0, nrow=nsim, ncol=length(initialPoint))
    regenProbsVec <- rep(0, times=nsim)
    bellVars <- rep(0, times=nsim)

    outMatrix[1,] <- initialPoint
    for(i in 2:nsim) {
      outMatrix[i,] <- genNextState(x=outMatrix[i-1,])
      if(regenProbsMissing) {
        regenProbsVec[i-1] <- smallFn(x=outMatrix[i-1,]) *
                              smallMeasure(y=outMatrix[i,]) /
                              transDens(x=outMatrix[i-1,], y=outMatrix[i,])
      }
      else {
        regenProbsVec[i-1] <- regenProbs(x=outMatrix[i-1,], y=outMatrix[i,])
      }
    }

    # generate regeneration points
    bellVars[-nsim] <- ifelse(runif(nsim-1) < regenProbsVec[-nsim], 1, 0)
    # cat("Done.\n")

    if(sum(bellVars) >= 2) {
      tourMatrix <- splitIntoTours(outMatrix, bellVars, numericParams)
      return(tourMatrix)
    }
    else {
      # cat("No complete tours found. Returning everything.\n")
      return(list(outMatrix, regenProbsVec))
    }
  }
  else { # Run for a specified number of tours
    # Look for first regeneration point
    # cat("Looking for first regeneration point .. ")
    findFirstPt <- 0
    while(findFirstPt == 0) {
      potentialNewY <- genNextState(x=initialPoint)
      if(regenProbsMissing) {
        regenProb <- smallFn(x=initialPoint) *
                     smallMeasure(y=potentialNewY) /
                     transDens(x=initialPoint, y=potentialNewY)
      }
      else {
        regenProb <- regenProbs(x=initialPoint, y=potentialNewY)
      }
      complete <- ifelse(runif(1) < regenProb, 1, 0)
      if(complete == 1) {
        findFirstPt <- 1
        newStart <- potentialNewY
      }
      else {
        initialPoint <- potentialNewY
      }
    }
    # cat("Done.\n")

    # Run the Markov chain for nTours tours
    # cat("Generating", nTours, "tours .. ", sep=" ")
    tourMatrix <- matrix(0, nrow=nTours, ncol=length(numericParams)+1)
    for(i in 1:nTours) {
      tmp <- genOneTour(newStart, numericParams)
      tourMatrix[i,] <- tmp$out
      newStart <- tmp$newStart
    }
    # cat("Done.\n")
    colnames(tourMatrix) <- c(paste("S", numericParams, sep=""), "N")

    return(tourMatrix)

  }

}

