# mpi.spawn.Rslaves(.......)

Y <- matrix(c(3.071, 3.871, 2.965,
              4.319, 4.396, 5.045,
              5.221, 4.876, 5.058,
              4.572, 5.116, 5.578, 
              5.351, 3.925, 4.217, 
              5.889, 4.893, 4.775, 
              5.192, 4.457, 5.097,
              4.477, 4.807, 5.345, 
              5.060, 5.271, 5.454, 
              5.188, 4.499, 5.340, 
              5.970, 5.660, 5.175, 
              5.619, 1.843, 5.545, 
              4.200, 5.294, 4.945), nrow=13, byrow=TRUE)
a <- -0.5; b <- 0; q <- 13

miVec <- rep(ncol(Y), q)
trtMeans <- rowMeans(Y)

# sigE, sigT, theta1, ... , thetaq, mu
initialPoint <- c(0.19,0.61,0.9*trtMeans,1)

Dmat <- matrix(c(0.1,0.2,0.8,0.9), nrow=2)

argList <- list(a=a, b=b, q=q, Y=Y, miVec=miVec, trtMeans = trtMeans)

initialRun <- runMarkovChainNoRegenS(genNextState = genNextStateOnewayDiffuse,
                                     genNextStateArgList = argList, 
                                     initialPoint=initialPoint, 
                                     nsim=33000)
distPt <- colMeans(initialRun)
regenProbsArgList <- list(xTil=distPt, Dmat=Dmat, a=a, b=b, Y=Y, q=q, miVec=miVec, 
                          trtMeans=trtMeans)

#tmp <- regenProbsFromRawOutput(initialRun, regenProbsOnewayDiffuse, regenProbsArgList)

out1 <- runMarkovChainRegenS(genNextState=genNextStateOnewayDiffuse,
                                 genNextStateArgList=argList,
                                 initialPoint = initialPoint, 
                                 regenProbs = regenProbsOnewayDiffuse,
                                 regenProbsArgList=regenProbsArgList,
                                 nTours=11298)

# These are for a parallel interactive run
# nTours <- 50551
#
#A <-  runMarkovChainRegenIP(genNextState = genNextStateOnewayDiffuse, 
#                      genNextStateArgList = argList, 
#                      initialPoint = initialPoint, 
#                      regenProbs = regenProbsOnewayDiffuse, 
#                      regenProbsArgList = regenProbsArgList, 
#                      nTours=nTours)
# mpi.close.Rslaves(dellog=FALSE)
# mpi.quit()
