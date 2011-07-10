# generate next state in  HLMM Gibbs sampler
genNextStateHLMM <- function(x, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n) {
  require(MASS)
  y <- rep(0, q+p+2)

  # generate lamD variable - recall that x is 
  # arranged as (lamD, lamR, u, beta)
  y[1] <- rgamma(1, q/2 + d1, rate=d2 + 0.5 * drop(crossprod(x[3:(2+q)])))

  # generate lamR
  y[2] <- rgamma(1, n/2 + r1, rate = r2 + 0.5 * 
                 drop(crossprod(Y - X %*% x[(2+q+1):(q+p+2)] - 
                 Z %*% x[3:(2+q)])))

  # generate xi
  Rmat <- y[2] * diag(n)
  Dmat <- y[1] * diag(q)
  tmpSigma1 <- cbind(t(Z) %*% Rmat %*% Z + Dmat, t(Z) %*% Rmat %*% X)
  tmpSigma2 <- cbind(t(X) %*% Rmat %*% Z, t(X) %*% Rmat %*% X + Bmat)
  covMat <- solve(rbind(tmpSigma1, tmpSigma2))

  meanVec <- drop(covMat %*% rbind(t(Z) %*% Rmat %*% Y,  
                              t(X) %*% Rmat %*% Y + Bmat %*% beta0))

  y[-(1:2)] <- mvrnorm(1, mu=meanVec, Sigma=covMat)

  y
}
