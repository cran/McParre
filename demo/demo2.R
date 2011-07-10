# True posterior means of \lam_theta, \theta_1 and \theta_2 are:
#                         0.1361      0.5208       -3.9258
# this is the data:
Y <- matrix(c(-1.87008415305704, -0.390234937174061, 0.284387730466376, 
              5.91202546805967, 2.08062321154896, -1.56760403993739, 
              -2.38159613905326, -4.63685644567788, -4.65962816209147, 
              -4.921806311828, -3.37217794128032, -6.04342597681985, 
              -2.94385203069467,-6.85034281162141), nrow=2, byrow=TRUE)

argList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, b1 = 4.3,
                a2 = 2.1, b2 = 20.3, Y=Y)

initialPoint <- rep(1, times=5)

xTil <- c(0.1363941,  0.1524292,  2.1420824,  0.5226216, -3.9227209)
Dmat <- matrix(0, nrow=3, ncol=2)
Dmat[1,] <- c(0.04446411, 0.2536926)
Dmat[2,] <- c(0.08605464, 0.2273910)
Dmat[3,] <- c(1.85707086, 2.4272405)

smallMeasureArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                            b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, 
                            xTil=xTil, D=Dmat)
smallFnArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                       b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, xTil=xTil, D=Dmat)

regenProbsArgList <- list(K = 2, M = 7, lam0 = 20.3, mu0 = 2.19, a1 = 2.1, 
                          b1 = 4.3, a2 = 2.1, b2 = 20.3, Y=Y, xTil=xTil, D=Dmat)

rm(Y, xTil, Dmat)

# Arguments have to be named!
runMarkovChainRegenS(genNextState = genNextStateOnewayPlain, 
                     genNextStateArgList = argList, 
                     initialPoint = initialPoint, 
                     transDens = transDensOnewayPlain, 
                     transDensArgList = argList, 
                     smallMeasure = smallMeasureOnewayPlain, 
                     smallMeasureArgList = smallMeasureArgList, 
                     smallFn = smallFnOnewayPlain, 
                     smallFnArgList = smallFnArgList, nTours=10) -> output1

computeCI(output1, 0.05)

coeffVar(output1)

runMarkovChainRegenS(genNextState = genNextStateOnewayPlain, 
                     genNextStateArgList = argList, 
                     initialPoint = initialPoint, 
                     regenProbs = regenProbsOnewayPlain, 
                     regenProbsArgList = regenProbsArgList, 
                     nTours=10) -> output1
