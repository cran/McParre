#include<R.h>
#include<Rmath.h>
#include<R_ext/BLAS.h>

void rdirichlet(double *shape1, int *dimension, double *output) {
  int i;
  double currentSum = 0.0;

  for(i=0; i<*dimension; i++) {
    output[i] = rgamma(shape1[i], 1.0);
    currentSum += output[i];
  }
  for(i=0; i<*dimension; i++) {
    output[i] = output[i]/currentSum;
  }
}

void rThetaGivenZ(int *zid, int *did, int *K, int *D, int *totalWords, 
                  double *alpha, double *output) {
  int i;
  double *outTable;

// set up a table for the counts and initialize to alpha
  outTable = (double *) Calloc((*D) * (*K), double);
  for(i=0; i<(*D)*(*K); i++)
    outTable[i] = alpha[i % (*K)];

// carry out the counting
  for(i=0; i<*totalWords; i++)
    outTable[ (did[i]-1)*(*K) + zid[i] - 1 ]++;
//  for(i=0; i<(*D)*(*K); i++)
//    Rprintf("%6.3f ", outTable[i]);

// generate theta values
  for(i=0; i<*D; i++)
    rdirichlet(outTable + i*(*K), K, output + i*(*K));

  Free(outTable);
}

void rBetaGivenZ(int *zid, int *wid, int *K, int *V, int *totalWords, 
                  double *betaPrior, double *output) {
  int i;
  double *outTable;

// set up a table for the counts and initialize to alpha
  outTable = (double *) Calloc((*V) * (*K), double);
  for(i=0; i<(*V)*(*K); i++)
    outTable[i] = betaPrior[i % (*V)];

//  carry out the counting
  for(i=0; i<*totalWords; i++)
    outTable[ (zid[i]-1)*(*V) + wid[i] - 1 ]++;
//  for(i=0; i<(*V)*(*K); i++)
//    Rprintf("%6.3f ", outTable[i]);
//  Rprintf("\n", outTable[i]);

// generate beta values
  for(i=0; i<*K; i++)
    rdirichlet(outTable + i*(*V), V, output + i*(*V));

  Free(outTable);
}

void rZidGivenBetaTheta(int *wid, int *did, int *K, int *D, int *V, 
                        int *totalWords, double *beta, double *theta, 
                        int *zOut) {
  int i, ii, *zOutTmp;
  double *normConst, *probVec;
  char noTrans = 'n';
  double alphaBLAS=1.0, betaBLAS=0.0;

  normConst = (double *) Calloc((*D)*(*V), double);
  probVec = (double *) Calloc(*K, double);
  zOutTmp = (int *) Calloc(*K, int);

  F77_CALL (dgemm) (&noTrans, &noTrans, V, D, K, &alphaBLAS, beta, V, 
                    theta, K, &betaBLAS, normConst, V);

//  GetRNGstate();
  for(i=0; i<*totalWords; i++){
    for(ii=0; ii<*K; ii++) {
      probVec[ii] = beta[ii*(*V) + wid[i]-1] * theta[(did[i]-1)*(*K) + ii] /
                    normConst[(did[i]-1)*(*V) + wid[i] - 1];
//      Rprintf("%8.3f ", probVec[ii]);
    }
    rmultinom(1, probVec, *K, zOutTmp);
    for(ii=0; ii<*K; ii++) {
      if(zOutTmp[ii] > 0) {
        zOut[i] = ii+1;
        break;
      }
    }
//    Rprintf("%02d\n", zOut[i]);
  }
//  PutRNGstate();

  Free(normConst);
  Free(probVec);
  Free(zOutTmp);
}

void genNextStateTopicModelInC (double *x, double *y, int *wid, int *did, 
                                 int *K, int *numDocs, int *V, int *totalWords, 
                                 double *alpha, double *betaPrior) {
  int i, *zidOld, *zidNew;
  
  zidOld = (int *) R_alloc(*totalWords, sizeof(int));
  zidNew = (int *) R_alloc(*totalWords, sizeof(int));
  for(i=0; i<*totalWords; i++)
    zidOld[i] = (int) x[i];
  
  GetRNGstate();

  // generate theta given z
  rThetaGivenZ(zidOld, did, K, numDocs, totalWords, alpha, y + *totalWords);

  // generate beta given z
  rBetaGivenZ(zidOld, wid, K, V, totalWords, betaPrior, 
              y + *totalWords + (*K)*(*numDocs));

  // generate z given beta and theta
  rZidGivenBetaTheta(wid, did, K, numDocs, V, totalWords, 
                     y + *totalWords + (*K)*(*numDocs), 
                     y + *totalWords, zidNew);
  for(i=0; i<*totalWords; i++)
    y[i] = (double) zidNew[i];

  PutRNGstate();

}
