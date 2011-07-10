# run McParre on MASS data
library(MASS)
library(McParre)

Y <- as.vector(petrol$Y)
n <- length(Y)

# Fixed hyperprior parameters
r1 <- 100
r2 <- 3.1

d1 <- 100
d2 <- 3

# Matrices 
X <- as.matrix(cbind(1, scale(petrol[,2:5], scale=F)))
Z <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=32, byrow=TRUE)

p <- ncol(X)
q <- ncol(Z)
Bmat <- 0.001*diag(p)
beta0 <- rep(0, p)

# initial point
x1 <- rep(1, p+q+2) 

argList <- list(Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, beta0=beta0, Bmat=Bmat, 
                p=p, q=q, n=n)

# run non-regenerative Markov Chain to get distinguished point and sets
runMarkovChainNoRegenS(genNextStateHLMM, argList, initialPoint=x1, 
                       nsim=10000) -> initialRun

initialRun <- initialRun[-(1:6000),]
xTil <- colMeans(initialRun)
newInitialPoint <- initialRun[4000,]

lamRInterval <- c(mean(initialRun[,2]) - 1.5*sd(initialRun[,2]), 
                  mean(initialRun[,2]) + 1.5*sd(initialRun[,2]))

lamDInterval <- c(mean(initialRun[,1]) - 1.5*sd(initialRun[,1]), 
                  mean(initialRun[,1]) + 1.5*sd(initialRun[,1]))

smallFnArgList <- list(xTil=xTil, Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, 
                       beta0=beta0, Bmat=Bmat, p=p, q=q, n=n, 
                       lamRInterval=lamRInterval,
                       lamDInterval=lamDInterval)

smallMeasureArgList <- list(xTil=xTil, Y=Y, X=X, Z=Z, r1=r1, r2=r2, d1=d2, 
                       beta0=beta0, Bmat=Bmat, p=p, q=q, n=n, 
                       lamRInterval=lamRInterval,
                       lamDInterval=lamDInterval)

runMarkovChainRegenS(genNextStateHLMM, argList, newInitialPoint, 
                     transDens=transDensHLMM, transDensArgList=argList, 
                     smallMeasure=smallMeasureHLMM, 
                     smallMeasureArgList=smallMeasureArgList, 
                     smallFn=smallFnHLMM, smallFnArgList=smallFnArgList, 
                     numericParams=1:(p+q+2), nTours=30) -> out1

# USE NAMED ARGUMENTS!
#runMarkovChainRegenP(genNextStateHLMM, argList, newInitialPoint, transDensHLMM, argList, smallMeasureHLMM, smallMeasureArgList, smallFnHLMM, smallFnArgList, numericParams=1:(p+q+2), nTours=30, prefix="rank", pathToLog="log/")
