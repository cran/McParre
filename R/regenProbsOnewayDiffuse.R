# regeneration probabilities for oneway model with diffuse priors.
# xTil is same length as x and y here, not like in Aixin's paper.
regenProbsOnewayDiffuse <- 
function(x, y, xTil, Dmat, a, b, Y, q, miVec, trtMeans) {

  # check if y is in distinguished set
  tmp <- apply(cbind(y[1:2], Dmat), 1, function(x)
               {findInterval(x[1], x[-1]) == 1})
  if(sum(tmp) < 2) return(0)

  # this part is to convert the sigma^2 parameter to lambda parametrization
  x[1:2] <- 1/x[1:2]
  y[1:2] <- 1/y[1:2]
  xTil[1:2] <- 1/xTil[1:2]

  one <- rep(1,q)
  M <- cbind(diag(q), -one)
  MM <- t(M) %*% M

  # remember, order is sigT, sigE, thetaVec, mu
  w1 <- t(x[-(1:2)]) %*% MM %*% x[-(1:2)] 
  w2 <- drop(crossprod(miVec, (trtMeans - x[3:(q+2)])^2))

  w1Star <- t(xTil[-(1:2)]) %*% MM %*% xTil[-(1:2)] 
  w2Star <- drop(crossprod(miVec, (trtMeans - xTil[3:(q+2)])^2))

  if(w1 > w1Star) 
    lat <- Dmat[1,2] else 
    lat <- Dmat[1,1]
  if(w2 > w2Star) 
    lae <- Dmat[2,2] else
    lae <- Dmat[2,1]

  logProb <- 0.5 * ( (w1 - w1Star)*(y[1] - lat) + 
                     (w2 - w2Star)*(y[2] - lae)  )

  exp(logProb)
}
