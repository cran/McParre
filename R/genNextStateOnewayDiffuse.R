# generate next state for oneway model with diffused priors.
genNextStateOnewayDiffuse <- 
function(x, a, b, Y, q, miVec, trtMeans) {

  # order is lamtheta, lame, theta, mu
  x[1:2] <- 1/x[1:2]

  one <- rep(1,q)
  M <- cbind(diag(q), -one)
  MM <- t(M) %*% M

  w10 <- t(x[-(1:2)]) %*% MM %*% x[-(1:2)] 
  w20 <- sum((Y - x[3:(q+2)])^2, na.rm=TRUE)

  lat <- rgamma(1, q/2 + a, rate=w10/2)
  lae <- rgamma(1, sum(miVec)/2 +b, rate=w20/2)

  d <- sqrt(lat + miVec*lae)
  t <- sum((lat*lae*miVec)/(lat + miVec*lae))
  invL <- diag(c(1/d, 1/sqrt(t)))
  invL[(q+1), 1:q] <- lat/(lat + miVec*lae)/sqrt(t)
  V <- t(invL) %*% invL

  zeta_zero <- V %*% c(miVec*trtMeans*lae, 0)
  newXi <- zeta_zero + t(invL) %*% rnorm(q+1)
 
  c(1/lat, 1/lae, drop(newXi))
}
