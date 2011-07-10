runMarkovChainNoRegenS <-
function(genNextState, genNextStateArgList=NULL, 
                                   initialPoint, nsim) {

# CHECKS ON genNextState AND ITS ARGUMENTS. ######################
  # check that the genNextState argument is a function.
  if(class(genNextState)!="function")
    stop("genNextState() has to be of 'function' class.\n")

  # check that genNextStateArgList is not missing and is a list.
  if(!is.null(genNextStateArgList)){
    if(class(genNextStateArgList)!="list")
      stop("genNextStateArgList has to be of 'list' class.\n")

    # After this, we should only have to call genNextState(x=...)
    formals(genNextState) <- c(genNextStateArgList, x=1)
  }

  # try running one step and stop if there are errors.
  cat("Checking to see if one-step generation function works .. ")
  firstPointGenerated <- genNextState(x=initialPoint)
  if(sum(is.na(firstPointGenerated)) > 0)
    stop("NA's generated in genNextState().\n")
  cat("OK\n")

  # check length of output argument
  if(length(firstPointGenerated)!=length(initialPoint))
    stop("Input/output states for genNextState() do not match in length.\n")
# END OF CHECKS ON genNextState AND ITS ARGUMENTS. ###############

# Run the Markov chain for nsim steps
  cat("Running Markov chain .. ")
  outMatrix <- matrix(0, nrow=nsim, ncol=length(initialPoint))

  outMatrix[1,] <- firstPointGenerated
  for(i in 2:nsim) {
    outMatrix[i,] <- genNextState(x=outMatrix[i-1,])
  }
  cat("Done.\n")

  outMatrix
}

