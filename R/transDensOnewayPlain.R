transDensOnewayPlain <- 
function(x, y, K, M, a1, b1, a2, b2, Y, lam0, mu0) {
  prod1 <- 1.0
  .C("transDensOnewayPlainInC", as.double(x), as.double(y), as.integer(K), 
     as.integer(M), as.double(a1), as.double(b1), as.double(a2), as.double(b2), 
     as.double(Y), as.double(lam0), as.double(mu0), as.double(prod1), 
     PACKAGE="McParre")[[12]]
}
