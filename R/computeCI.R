computeCI <-
function(outObj, alpha=0.05, paramsOfInterest) {
  # only compute the CI for the parameters of interest.
  if(!missing(paramsOfInterest)) {
    if(max(paramsOfInterest) > ncol(outObj)-1)
      stop("parameter indices incorrectly specified.\n")
    outObj <- outObj[,c(paramsOfInterest,ncol(outObj))]
  }

  NcolPos <- ncol(outObj)
  R <- nrow(outObj)  # total number of tours
  sumNi <- sum(outObj[ , NcolPos])  # total number of iterations over all tours.
  gBar <- colSums(outObj)[-NcolPos]/sumNi  # ergodic average

  sigMat <- apply(outObj, 1, function(xi) {(xi[-NcolPos] - xi[NcolPos]*gBar)^2})
  
  sigVec <- (rowSums(sigMat)*R)/(sumNi^2)

  outMat <- matrix(0, nrow=ncol(outObj)-1, ncol=3)
  rownames(outMat) <- colnames(outObj)[-NcolPos]
  colnames(outMat) <- c('lower', 'ptEstimate', 'upper')
  outMat[,1] <- gBar - qnorm(alpha/2, lower.tail=FALSE)*sqrt(sigVec/R)
  outMat[,2] <- gBar
  outMat[,3] <- gBar + qnorm(alpha/2, lower.tail=FALSE)*sqrt(sigVec/R)

  outMat

}

