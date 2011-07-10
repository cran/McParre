regenProbsFromRawOutput <-
function(outMatrix, regenProbs, regenProbsArgList=NULL, 
         transDens, transDensArgList=NULL,
         smallMeasure, smallMeasureArgList=NULL,
         smallFn, smallFnArgList=NULL) {

  # Check that either regenProbs, or all of smallFn, smallMeasure and 
  # transDens are supplied.
  regenProbsMissing <- missing(regenProbs)
  if(regenProbsMissing){
    transDensMissing <- missing(transDens)
    smallFnMissing <- missing(smallFn)
    smallMeasureMissing <- missing(smallMeasure)
    if(transDensMissing || smallFnMissing || smallMeasureMissing)
      stop("Either regenProbs(), or ALL of transDens(), smallFn() and 
            smallMeasure() have to be supplied.\n")
  }

  # Check the supplied functions and arguments.

  if(regenProbsMissing) {
# CHECKS ON transDens AND ITS ARGUMENTS. ######################
  # check that the genNextState argument is a function.
  if(class(transDens)!="function")
    stop("transDens() has to be of 'function' class.\n")

  # check that transDensArgList is not missing and is a list.
  if(!is.null(transDensArgList)){
    if(class(transDensArgList)!="list")
      stop("transDensArgList has to be of 'list' class.\n")

    # After this, we should only have to call genNextState(x=...)
    formals(transDens) <- c(transDensArgList, x=1, y=1)
  }

# CHECKS ON smallMeasure AND ITS ARGUMENTS. ######################
  # check that the smallMeasure argument is a function.
  if(class(smallMeasure)!="function")
    stop("smallMeasure() has to be of 'function' class.\n")

  # check that smallMeasureArgList is not missing and is a list.
  if(!is.null(smallMeasureArgList)){
    if(class(smallMeasureArgList)!="list")
      stop("smallMeasureArgList has to be of 'list' class.\n")

    # After this, we should only have to call genNextState(x=...)
    formals(smallMeasure) <- c(smallMeasureArgList, y=1)
  }

# CHECKS ON smallFn AND ITS ARGUMENTS. ######################
  # check that the smallFn argument is a function.
  if(class(smallFn)!="function")
    stop("smallFn() has to be of 'function' class.\n")

  # check that smallFnArgList is not missing and is a list.
  if(!is.null(smallFnArgList)){
    if(class(smallFnArgList)!="list")
      stop("smallFnArgList has to be of 'list' class.\n")

    # After this, we should only have to call genNextState(x=...)
    formals(smallFn) <- c(smallFnArgList, x=1)
  }
  }
  else {
# CHECKS ON regenProbs AND ITS ARGUMENTS. ######################
  # check that the regenProbs argument is a function.
  if(class(regenProbs)!="function")
    stop("regenProbs() has to be of 'function' class.\n")

  # When regenProbsArgList is not NULL, check that it is a list.
  if(!is.null(regenProbsArgList)){
    if(class(regenProbsArgList)!="list")
      stop("regenProbsArgList has to be of 'list' class.\n")

    # After this, we should only have to call genNextState(x=...)
    formals(regenProbs) <- c(regenProbsArgList, x=1, y=1)
  }
  }

  nsim <- nrow(outMatrix)
  regenProbsVec <- rep(NA, times=nsim)
  for(i in 2:nsim) {
    if(regenProbsMissing) {
      regenProbsVec[i-1] <- smallFn(x=outMatrix[i-1,]) * 
                            smallMeasure(y=outMatrix[i,]) /
                            transDens(x=outMatrix[i-1,], y=outMatrix[i,])
    }
    else {
      regenProbsVec[i-1] <- regenProbs(x=outMatrix[i-1,], y=outMatrix[i,])
    }
  }

  regenProbsVec

}

