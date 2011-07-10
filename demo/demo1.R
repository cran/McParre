# true values:    0.1357   -3.8330    0.1888
Y <- matrix(c(-4.60454421431755,-4.87571538076115,-8.45342137358092,
              -5.7527706056414,-3.31470408208705,-0.755164672373152,
              -4.73149297074921,1.29639133204756,0.095231285316642,
              -3.50505074286717,-0.88748529070148,-0.831152640489548,
              0.0802442700977176,3.23454790727842), nrow=2, byrow=TRUE)
argList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, b1 = 4.3,
                a2 = 2.1, b2 = 20.3, Y=Y)

initialPoint <- rep(1, times=5)

initialRun <- runMarkovChainNoRegenS(genNextState = genNextStateOnewayPlain,
                                     genNextStateArgList = argList, 
                                     initialPoint=initialPoint, 
                                     nsim=5000)

xTil <- colMeans(initialRun)
Dmat <- matrix(0, nrow=3, ncol=2)
Dmat[1,] <- c(quantile(initialRun[,1], 0.1), quantile(initialRun[,1], 0.9))
Dmat[2,] <- c(quantile(initialRun[,2], 0.1), quantile(initialRun[,2], 0.9))
Dmat[3,] <- c(quantile(initialRun[,3], 0.1), quantile(initialRun[,3], 0.9))

smallMeasureArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                            b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, 
                            xTil=xTil, D=Dmat)

smallFnArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                       b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, xTil=xTil, D=Dmat)

regenProbsArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                          b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, xTil=xTil, D=Dmat)

rm(Y, xTil, Dmat)

regenerationProbs <- regenProbsFromRawOutput(outMatrix=initialRun, 
                      transDens = transDensOnewayPlain, 
                      transDensArgList = argList, 
                      smallMeasure = smallMeasureOnewayPlain, 
                      smallMeasureArgList = smallMeasureArgList, 
                      smallFn = smallFnOnewayPlain, 
                      smallFnArgList = smallFnArgList)

regenerationProbs2 <- regenProbsFromRawOutput(outMatrix=initialRun, 
                       regenProbs = regenProbsOnewayPlain,
                       regenProbsArgList = regenProbsArgList)
