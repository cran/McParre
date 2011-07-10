# compute the transition density. note that it is not the full density
# we only need the first lam part (see equation (15) of scand journal paper)
smallMeasureHLMM <- function(xTil, y, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, 
                             q, n, lamRInterval, lamDInterval) {
  tmp1 <- findInterval(y[1], lamDInterval)
  tmp2 <- findInterval(y[2], lamRInterval)
  if((tmp1!=1)||(tmp2!=1))
    return(0)

  dgamma(y[1], q/2 + d1, rate = d2 + 0.5 * drop(crossprod(xTil[3:(2+q)]))) * 
  dgamma(y[2], n/2 + r1, rate = r2 + 0.5 * 
                 drop(crossprod(Y - X %*% xTil[(2+q+1):(q+p+2)] - 
                 Z %*% xTil[3:(2+q)])))
}
