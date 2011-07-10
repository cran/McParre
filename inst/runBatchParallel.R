# mpi.spawn.Rslaves(.......)
library(MASS)
library(McParre, lib.loc="/scratch/ufhpc/viknesh/Rlibs")

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

distPt <- c(0.1855264, 0.6197913, 4.2301059, 4.7243763, 4.9007954, 4.9158586, 
            4.6890765, 4.9543741, 4.8508905, 4.8341545, 4.9821088, 4.8858458, 
            5.1128051, 4.6234750, 4.8122107, 4.8094787)

regenProbsArgList <- list(xTil=distPt, Dmat=Dmat, a=a, b=b, Y=Y, q=q, 
                          miVec=miVec, trtMeans=trtMeans)

runMarkovChainRegenP(genNextState=genNextStateOnewayDiffuse,
                     genNextStateArgList=argList,
                     initialPoint = initialPoint, 
                     regenProbs = regenProbsOnewayDiffuse,
                     regenProbsArgList=regenProbsArgList, 
                     prefix="rank", pathToLog="log/", nTours=1120)
 mpi.quit()

