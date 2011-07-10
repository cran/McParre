smallMeasureOnewayPlain <- 
function(xTil, y, D, K, M, a1, b1, a2, b2, Y, lam0, mu0) {

  tmp <- apply(cbind(y[1:3], D), 1, function(x) 
               {findInterval(x[1], x[-1]) == 1})
  if(sum(tmp) < 3) return(0)

  prod1 <- 1.0
  .C("transDensOnewayPlainInC", as.double(xTil), as.double(y), as.integer(K), 
     as.integer(M), as.double(a1), as.double(b1), as.double(a2), as.double(b2), 
     as.double(Y), as.double(lam0), as.double(mu0), as.double(prod1), 
     PACKAGE="McParre")[[12]]
}
