#include<R.h>
#include<Rmath.h>

void genNextStateOnewayPlainInC(double *x, int *K, int *M, double *a1, 
                                double *b1, double *a2, double *b2, 
                                double *Y, double *lam0, double *mu0,
                                double *y) {
  double sum1=0.0, sum2=0.0, sd1=0.0, meanTmp=0.0, *meanY;
  int i, ii;

  meanY = (double *) R_alloc(*K, sizeof(double));
  for(ii=0; ii<*K; ii++)
    meanY[ii] = 0.0;

  GetRNGstate();

  for(ii=0; ii<*K; ii++) {
    sum1 += pow(x[3+ii] - x[2], 2) ;
    meanTmp += x[3+ii] * 1.0/(*K);

    for(i=0; i<(*M); i++) {
      sum2 += pow(Y[i*(*K)+ii] - x[3+ii], 2) ;
      meanY[ii] += Y[i*(*K)+ii] * 1.0/(*M);
    }
  }

  y[0] = rgamma(*K/2.0 + *a1, 1.0/(*b1 + 0.5 * sum1));

  y[1] = rgamma((*K)*(*M)/2.0 + *a2, 1.0/(*b2 + 0.5 * sum2));

  y[2] = rnorm((*lam0 * (*mu0) + (*K)*y[0]*meanTmp)/(*lam0 + *K * y[0]),
               sqrt(1.0/(*lam0 + (*K)*y[0])));
 
  sd1 = sqrt(1.0/(y[0] + *M * y[1]));
  for(ii=3; ii<*K+3; ii++)
    y[ii] = rnorm((y[0]*y[2] + *M * y[1] * meanY[ii-3])/(y[0] + *M * y[1]),
                  sd1);

  PutRNGstate();

}

void transDensOnewayPlainInC(double *x, double *y, int *K, int *M, double *a1, 
                             double *b1, double *a2, double *b2, 
                             double *Y, double *lam0, double *mu0, 
                             double *prod1) {

  double sum1=0.0, sum2=0.0, sd1=0.0, meanTmp=0.0, *meanY;
  int i, ii;

  meanY = (double *) R_alloc(*K, sizeof(double));
  for(ii=0; ii<*K; ii++)
    meanY[ii] = 0.0;

  for(ii=0; ii<*K; ii++) {
    sum1 += pow(x[3+ii] - x[2], 2) ;
    meanTmp += x[3+ii] * 1.0/(*K);

    for(i=0; i<(*M); i++) {
      sum2 += pow(Y[i*(*K)+ii] - x[3+ii], 2) ;
      meanY[ii] += Y[i*(*K)+ii] * 1.0/(*M);
    }
  }

  *prod1 = dgamma(y[0], *K/2.0 + *a1, 1.0/(*b1 + 0.5 * sum1), 0) *
          dgamma(y[1], (*K)*(*M)/2.0 + *a2, 1.0/(*b2 + 0.5 * sum2), 0) *
          dnorm(y[2], (*lam0 * (*mu0) + (*K)*y[0]*meanTmp)/(*lam0 + *K * y[0]),
                sqrt(1.0/(*lam0 + (*K)*y[0])), 0);
 
  sd1 = sqrt(1.0/(y[0] + *M * y[1]));
  for(ii=3; ii<*K+3; ii++)
    *prod1 *= dnorm(y[ii], (y[0]*y[2] + *M * y[1] * meanY[ii-3])/
                    (y[0] + *M * y[1]), sd1, 0);

}
