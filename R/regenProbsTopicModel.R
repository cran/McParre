regenProbsTopicModel <- 
function(x, y, xTil, distSetThetaAll, distSetBetaAll, K, V, 
         wid, did, totalWords, numDocs)
{
# check if points are in distinguished set. Otherwise return 0
  for(ii in 1:numDocs) {
    theta1 <- y[(totalWords + 1 + (ii-1)*K):(totalWords + ii*K)]
    #tmpVert <- verticesFromDistSet(distSetThetaAll[(1+(ii-1)*K):(ii*K),])
    tmp <- checkIfInDistSet(theta1, distSetThetaAll[(1+(ii-1)*K):(ii*K),])
    if(!tmp) return(0)
  }
  for(ii in 1:K) {
    beta1 <- y[(totalWords + 1 + numDocs*K + (ii-1)*V):
               (totalWords + numDocs*K + ii*V)]
    tmp <- checkIfInDistSet(beta1, distSetBetaAll[(1+(ii-1)*V):(ii*V),])
    if(!tmp) return(0)
  }

  smallFn <- 0
  zidX <- x[1:totalWords]
  zidTilde <- xTil[1:totalWords]
  njZidX <- as.vector(unlist(tapply(zidX, did, tabulate, nbins=K)))
  njZidTilde <- as.vector(unlist(tapply(zidTilde, did, tabulate, nbins=K)))
  njDiff <- ifelse(njZidX - njZidTilde > 0, njZidX - njZidTilde, 0)
  for(j in 1:numDocs) {
    coeff1 <- njDiff[((j-1)*K + 1):(j*K)]
    smallFn <- smallFn + smallFnComponent(coeff1, 
                verticesFromDistSet(distSetThetaAll[(1+(j-1)*K):(j*K),]))
  }
 
  MjiZidX <- matrix(unlist(tapply(zidX, wid, tabulate, nbins=K)), nrow=K)
  MjiZidTilde <- matrix(unlist(tapply(zidTilde, wid, tabulate, nbins=K)), 
                        nrow=K)
  MjiDiff <- ifelse(MjiZidX - MjiZidTilde > 0, MjiZidX - MjiZidTilde, 0)
  for(j in 1:K) {
    coeff1 <- MjiDiff[j,]
    smallFn <- smallFn + smallFnComponent(coeff1, 
                verticesFromDistSet(distSetBetaAll[(1+(j-1)*V):(j*V),]))
  }
  betaY <- matrix(y[(1+totalWords+K*numDocs):length(y)], 
                  nrow=K, ncol=V, byrow=TRUE)
  A1 <- sum((MjiZidTilde  - MjiZidX) * log(betaY)) # log of the p2a ratio

  thetaY <- matrix(y[(totalWords+1):(totalWords+K*numDocs)], nrow=K)
  A2 <- sum((njZidTilde  - njZidX) * log(thetaY)) # log of the p2a ratio

  exp(smallFn + A1 + A2)
}
