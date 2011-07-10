splitIntoTours <-
function(outMatrix, bellVars, numericParams) {
  numCompleteTours <- sum(bellVars) - 1
  if(numCompleteTours <= 0) return(NA)

  tourEndPts <- which(bellVars == 1)  
  tourLengths <- diff(tourEndPts)     # lengths of completed tours
  tourIdentifiers <- rep(1:numCompleteTours, times=tourLengths)

  subMatrix <- outMatrix[(tourEndPts[1]+1):tourEndPts[numCompleteTours+1], 
                         numericParams]
  cat("Splitting raw output into tours .. ")
  l <- split(as.data.frame(subMatrix), tourIdentifiers)
  tourMatrix <- unlist(lapply(l, colSums))
  tourMatrix <- matrix(tourMatrix, ncol=ncol(subMatrix), byrow=TRUE)
  tourMatrix <- cbind(tourMatrix, tourLengths)
  colnames(tourMatrix) <- c(paste("S", numericParams, sep=""), "N")
  cat("Done.\n")
  tourMatrix
}

