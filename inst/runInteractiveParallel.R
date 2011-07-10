# spawn slaves with mpi.spawn.Rslaves, or using provided RprofileINT
library(McParre)
mpi.remote.exec(library(McParre, lib.loc="/scratch/ufhpc/viknesh/Rlibs"))

demo(demo2)
nsim <- 150
nTours <- 35

A <-  runMarkovChainRegenIP(genNextState = genNextStateOnewayPlain, 
                      genNextStateArgList = argList, 
                      initialPoint = initialPoint, 
                      regenProbs = regenProbsOnewayPlain, 
                      regenProbsArgList = regenProbsArgList, 
                      nTours=nTours)

#A <-  runMarkovChainRegenIP(genNextState = genNextStateOnewayPlain, 
#                      genNextStateArgList = argList, 
#                      initialPoint = initialPoint, 
#                      regenProbs = regenProbsOnewayPlain, 
#                      regenProbsArgList = regenProbsArgList, 
#                      nsim=nsim)

# mpi.close.Rslaves(dellog=FALSE)
# mpi.quit()
