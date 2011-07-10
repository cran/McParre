genNextStateOnewayPlain <- 
function(x, K, M, a1, b1, a2, b2, Y, lam0, mu0) {
  y <- rep(0, times=K+3)

  .C("genNextStateOnewayPlainInC", as.double(x), as.integer(K), as.integer(M),
     as.double(a1), as.double(b1), as.double(a2), as.double(b2), 
     as.double(Y), as.double(lam0), as.double(mu0), as.double(y), 
     PACKAGE="McParre")[[11]]
}
