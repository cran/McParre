# compute the small function
smallFnHLMM <- function(x, xTil, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, 
                        n, lamRInterval, lamDInterval) {
  
  ux <- x[3:(2+q)]
  uxTil <- xTil[3:(2+q)]

  fac1 <- drop(crossprod(ux) - crossprod(uxTil))
  if (fac1 > 0)
    rat1 <- lamDInterval[2] else 
    rat1 <- lamDInterval[1]
  tmp1 <- (d2 + 0.5 * drop(crossprod(ux)))/(d2 + 0.5 * drop(crossprod(uxTil)))
  tmp1 <- tmp1^(q/2 + d1) * exp(-0.5*fac1*rat1)

  vx <- Y - X %*% x[(2+q+1):(q+p+2)] - Z %*% x[3:(2+q)]
  vxTil <- Y - X %*% xTil[(2+q+1):(q+p+2)] - Z %*% xTil[3:(2+q)]

  fac2 <- drop(crossprod(vx) - crossprod(vxTil))
  if (fac2 > 0)
    rat2 <- lamRInterval[2] else 
    rat2 <- lamRInterval[1]
  tmp2 <- (r2 + 0.5 * drop(crossprod(vx)))/(r2 + 0.5 * drop(crossprod(vxTil)))
  tmp2 <- tmp2^(n/2 + r1) * exp(-0.5*fac2*rat2)

  tmp1 * tmp2
}
