smallFnOnewayPlain <-
function(x, xTil, D, K, M, a1, b1, a2, b2, Y, lam0, mu0) {
  d1 <- D[1,]
  d2 <- D[2,]
  d3 <- D[3,]

  findWhichCase <- function() {
    if ((rt1<d1[1])&(rt2<d1[1]))
      return(1)
    else if ((rt1<d1[1])&(rt2<d1[2]))
      return(2)
    else if ((rt1<d1[1])&(rt2>d1[2]))
      return(3)
    else if ((rt1>d1[1])&(rt2<d1[2]))
      return(4)
    else if ((rt1<d1[2])&(rt2>d1[2]))
      return(5)
    else
      return(6)
  }
  s1 <- function(z1, z2, z3) {
    v1 <- 1/(lam0 + K*z1)
    m1.top <- (lam0*mu0 + K*z1*mean(x[-( 1:3 )]))/(lam0 + K*z1)
    m1.bot <- (lam0*mu0 + K*z1*mean(xTil[-( 1:3 )]))/(lam0 + K*z1)

    p1a <- dgamma(z1, shape = A1, rate = B1.top)
    p2a <- dgamma(z2, shape = A2, rate = B2.top )
    p3a <- dnorm(z3, mean = m1.top, sd =  sqrt(v1))
    p1b <- dgamma(z1, shape = A1, rate = B1.bot)
    p2b <- dgamma(z2, shape = A2, rate = B2.bot )
    p3b <- dnorm(z3, mean = m1.bot, sd =  sqrt(v1))
    (p1a * p2a * p3a)/(p1b * p2b * p3b)
  }

  B1.top <- b1 + 0.5 * sum( (x[-( 1:3 )] - x[3])^2 )
  B1.bot <- b1 + 0.5 * sum( (xTil[-( 1:3 )] - xTil[3])^2 )
  C1 <-  B1.top - B1.bot
  A1 <- ( K/2 + a1 )
    
  # Set z[2]
  B2.top <- b2 + 0.5 * sum( (Y - x[-( 1:3 )])^2 )
  B2.bot <- b2 + 0.5 * sum( (Y - xTil[-( 1:3 )])^2 )
  C2 <-  B2.top - B2.bot
  A2 <- (0.5)*M*K + a2
  if (C2 < 0) 
    z2 <- d2[1]
  else
    z2 <- d2[2]

  # Set z[3]
  m2.top <- mean(x[-(1:3)])
  m2.bot <- mean(xTil[-(1:3)])
  if (m2.bot - m2.top > 0) 
    z3 <- d3[2]
  else 
    z3 <- d3[1]
  
  coeff01 <- (0.5)*(K^2)*(-2*C1 - K*(m2.bot-m2.top)*(2*z3 - m2.bot - m2.top)) 
  coeff02 <- 2*lam0*K*(-C1 - K*z3*(m2.bot - m2.top))
             + K*(m2.bot - m2.top)*(lam0*K*(m2.top + m2.bot))
  coeff03 <- (lam0^2)*(-C1 - K*z3*(m2.bot - m2.top) + K*mu0*(m2.bot-m2.top))
  
  discriminant <- coeff02^2 - 4*coeff01*coeff03
  if (discriminant > 0) {
    (-coeff02 + sqrt(coeff02^2 - 4*coeff01*coeff03))/(2*coeff01) -> rt2
    (-coeff02 - sqrt(coeff02^2 - 4*coeff01*coeff03))/(2*coeff01) -> rt1
    if (rt2 < rt1) {
      tmp <- rt1
      rt1 <- rt2
      rt2 <- tmp
    }
  
    ret <- switch(findWhichCase(), 
           min(s1(d1[1],z2,z3), s1(d1[2],z2,z3)),
           min(s1(d1[1],z2,z3), s1(rt2, z2,z3), s1(d1[2],z2,z3)),
           min(s1(d1[1],z2,z3), s1(d1[2],z2,z3)),
           min(s1(d1[1],z2,z3), s1(rt1, z2,z3), s1(rt2, z2,z3), 
               s1(d1[2],z2,z3)),
           min(s1(d1[1],z2,z3), s1(rt1, z2,z3), s1(d1[2],z2,z3)),
           min(s1(d1[1],z2,z3), s1(d1[2],z2,z3)))
  }
  else 
    ret <- min(s1(d1[1],z2,z3), s1(d1[2],z2,z3))
  ret
}

