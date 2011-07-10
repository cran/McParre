# Generate DATA
 set.seed(1)
 V <- 3
 numDocs <- 2
 K <- 2
 nD <- 10 # number of words per document
 
 betaMat <- matrix(c(0.80, 0.10, 0.10,
                     0.07, 0.43, 0.53), 
                   nrow=K, byrow=TRUE)
 if(nrow(betaMat) != K)
   stop("check nrow(betaMat)")
 if(ncol(betaMat) != V)
   stop("check ncol(betaMat)")

 thetaMat <- matrix(c(0.05, 0.95,
                      0.95, 0.05), nrow=K)
 if(nrow(thetaMat) != K)
   stop("check nrow(thetaMat)")
 if(ncol(thetaMat) != numDocs)
   stop("check ncol(thetaMat)")

 wid <- rep(0, numDocs * nD)
 did <- rep(1:numDocs, each=nD)

 for(i in 1:length(wid)) {
   pickTopic <- sample(K, size=1, prob=thetaMat[ , did[i]])
   pickWord <- sample(V, size=1, prob=betaMat[pickTopic, ])
   wid[i] <- pickWord
 }
# Data generated

# Run Gibbs sampler for initial Run
 nsim <- 50000

 zid <- sample(K, length(did), replace=TRUE)
 betaMat <- matrix(abs(rnorm(K*V)), nrow=K, ncol=V)
 betaMat <- betaMat/rowSums(betaMat)

 thetaMat <- matrix(abs(rnorm(numDocs*K)), nrow=K)
 thetaMat <- thetaMat/rep(colSums(thetaMat), each=K)

 alpha <- rep(1.1, times=K)
 betaPrior <- rep(1.3, times=V)

 totalWords <- length(wid)

# generate points
 xState <- c(zid, as.vector(thetaMat), as.vector(t(betaMat)))
 argList <- list(wid=wid, did=did, K=K, numDocs=numDocs, V=V, 
                 totalWords=totalWords, alpha=alpha, betaPrior=betaPrior)
 
 outMat <- runMarkovChainNoRegenS(genNextStateTopicModel, argList, 
                                  initialPoint=xState, nsim=nsim)
# Finished running Gibbs sampler

# Get distinguished sets and points
constructDistSet03 <- function(dataMx, prob) {
  d <- ncol(dataMx)
  Amat <- matrix(0, nrow=d, ncol=d-1)
  dVec <- rep(0, times=d)

  for(ii in 1:(d-1)) {
    Amat[ii, ii] <- -1
    dVec[ii] <- -quantile(dataMx[,ii], prob)
  }
  Amat[d,] <- 1
  dVec[d] <- quantile(rowSums(dataMx[ , 1:(d-1)]), 1-prob)
  if( sum(dVec[-d]) >= dVec[d] )
    stop("Error constructing distinguished set\n")

  return(cbind(Amat, dVec))
}

 beta1RV <- outMat[,(1+totalWords+K*numDocs):(V+totalWords+K*numDocs)]
 beta2RV <- outMat[,(1+V+totalWords+K*numDocs):(2*V+totalWords+K*numDocs)]
 distSetBeta1 <- constructDistSet03(beta1RV, 0.3)
 distSetBeta2 <- constructDistSet03(beta2RV, 0.3)
 distSetBetaAll <- rbind(distSetBeta1, distSetBeta2)

 thetaRV <- outMat[,(totalWords+1):(totalWords+K*numDocs)]
 thetaTilde <- colMeans(thetaRV)
 betaTilde1 <- colMeans(beta1RV)
 betaTilde2 <- colMeans(beta2RV)
 zidTilde <- round(colMeans(outMat[,1:totalWords]))
 xTil <- c(zidTilde, thetaTilde, betaTilde1, betaTilde2)

 thetaDistSet <- matrix(0, nrow=numDocs, ncol=2)
 for(i in 1:numDocs) {
   tmp <- outMat[ , totalWords + (i-1)*K + 1]
   thetaDistSet[i,] <- c(exp(mean(log(tmp)) - 0.7*sd(log(tmp))),
                         exp(mean(log(tmp)) + 0.7*sd(log(tmp))))
 }
 thetaDistSet <- cbind(-1*thetaDistSet[,1], thetaDistSet[,2])
 distSetThetaAll <- cbind(1, as.vector(t(thetaDistSet)))
 distSetThetaAll[c(1,3,5),1] <- -1
 rm(tmp, thetaDistSet, thetaRV, zidTilde, betaTilde2, betaTilde1, thetaTilde)
 rm(beta1RV, beta2RV, distSetBeta1, distSetBeta2)
# distinguished sets/points found

# run regeneration
 regenProbsArgList <- list(xTil=xTil, distSetThetaAll=distSetThetaAll, 
                           distSetBetaAll=distSetBetaAll, K=K, V=V, wid=wid, 
                           did=did, totalWords=totalWords, numDocs=numDocs)

# tmp <- regenProbsFromRawOutput(outMat, regenProbsTopicModel, 
#                                regenProbsArgList = regenProbsArgList)

 newInitialPoint <- drop(outMat[nsim,])
 runMarkovChainRegenS(genNextState=genNextStateTopicModel, 
                      genNextStateArgList=argList,
                      initialPoint = newInitialPoint,
                      regenProbs = regenProbsTopicModel,
                      regenProbsArgList = regenProbsArgList,
                      numericParams = (totalWords+1):length(newInitialPoint),
                      nTours=2) -> outMat2
