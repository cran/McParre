n <- 6
# Fixed hyperprior parameters
r1 <- 3.2
r2 <- 3.1

d1 <- 1.2
d2 <- 2.1

p <- 3
q <- n
Bmat <- diag(p)

beta0 <- 1:p

# Matrices
Z <- diag(n)
X <- matrix(0, nrow=n, ncol=p)
X[1:2,1] <- 1
X[3:4,2] <- 1
X[5:6,3] <- 1

Y <- c(3.033603,  2.506486,  2.883229,  5.675019,  1.487337, -0.831464)

argList <- list(Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, beta0=beta0, Bmat=Bmat, 
                p=p, q=q, n=n)

## run non-regenerative Markov Chain to get distinguished point and sets
#runMarkovChainNoRegenS(genNextState = genNextStateHLMM, 
#                       genNextStateArgList = argList, initialPoint=x1, 
#                       nsim=10000) -> initialRun

newInitialPoint <- c(1.0146049, 0.7927502, 1.3970720, -0.3445606, -1.2445079, 
                     1.2942851, -0.5648055, -1.5207528, 1.6624752, 3.5163871, 
                     2.2278038)

xTil <- c(0.867455419, 0.934927379, 0.719274587, 0.406985243, -0.001597925,
          1.470935324, -0.255071684, -1.449298429, 1.719086371, 2.931709227,
          1.916930803)

# set up distinguished sets
lamDInterval <- c(-0.03228341,  1.76719425) 
lamRInterval <- c(0.1344841, 1.7353707)

smallFnArgList <- list(xTil=xTil, Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, 
                       beta0=beta0, Bmat=Bmat, p=p, q=q, n=n, 
                       lamRInterval=lamRInterval,
                       lamDInterval=lamDInterval)

smallMeasureArgList <- list(xTil=xTil, Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, 
                       beta0=beta0, Bmat=Bmat, p=p, q=q, n=n, 
                       lamRInterval=lamRInterval,
                       lamDInterval=lamDInterval)

runMarkovChainRegenS(genNextState = genNextStateHLMM, 
                     genNextStateArgList = argList, 
                     initialPoint = newInitialPoint, 
                     transDens = transDensHLMM, 
                     transDensArgList = argList, 
                     smallMeasure = smallMeasureHLMM, 
                     smallMeasureArgList = smallMeasureArgList, 
                     smallFn = smallFnHLMM, 
                     smallFnArgList = smallFnArgList, 
                     numericParams=1:(p+q+2), nsim=1000) -> output2

#runMarkovChainRegenS(genNextState = genNextStateHLMM, 
#                     genNextStateArgList = argList, 
#                     initialPoint = newInitialPoint, 
#                     transDens = transDensHLMM, 
#                     transDensArgList = argList, 
#                     smallMeasure = smallMeasureHLMM, 
#                     smallMeasureArgList = smallMeasureArgList, 
#                     smallFn = smallFnHLMM, 
#                     smallFnArgList = smallFnArgList, 
#                     numericParams=1:(p+q+2), nTours=1000) -> output2
#
# For parallel BATCH Mode - remember to use the correct .Rprofile file!
#
#runMarkovChainRegenP(genNextState = genNextStateHLMM, 
#                     genNextStateArgList = argList, 
#                     initialPoint = newInitialPoint, 
#                     transDens = transDensHLMM, 
#                     transDensArgList = argList, 
#                     smallMeasure = smallMeasureHLMM, 
#                     smallMeasureArgList = smallMeasureArgList, 
#                     smallFn = smallFnHLMM, 
#                     smallFnArgList = smallFnArgList, 
#                     numericParams=1:(p+q+2), nTours=43, 
#                     prefix="rank", pathToLog="log/")
#mpi.quit()
