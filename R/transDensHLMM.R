# compute the transition density. note that it is not the full density
# we only need the first lam part (see equation (15) of scand journal paper)
transDensHLMM <- function(x, y, Y, X, Z, r1, r2, d1, d2, beta0, Bmat, p, q, n) {

  dgamma(y[1], q/2 + d1, rate = d2 + 0.5 * drop(crossprod(x[3:(2+q)]))) * 
  dgamma(y[2], n/2 + r1, rate = r2 + 0.5 * 
                 drop(crossprod(Y - X %*% x[(2+q+1):(q+p+2)] - 
                 Z %*% x[3:(2+q)])))
}
