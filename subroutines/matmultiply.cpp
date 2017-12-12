
#include "matmultiply.h"


void matmultiply_nomkl(int nmol, int ncomplex, double *indivconc, double *complexconc, double *A)
{
  int j, i;
  for(j=0;j<nmol;j++)
    indivconc[j]=0;
  
  
  for(j=0;j<nmol;j++){
    for(i=0;i<ncomplex;i++){
      indivconc[j]+=A[i*nmol+j]*complexconc[i];
    }
  }
  
}


void matmultiply(int nmol, int ncomplex, double *indivconc, double *complexconc, double *A)
{
  int j, i;
  for(j=0;j<nmol;j++)
    indivconc[j]=0;
  
  char trans='N';
  int M=nmol;
  int N=ncomplex;
  double alpha=1;
  int lda=M;
  int incx=1;
  double beta=0;
  
  dgemv(&trans, &M, &N, &alpha, A, &lda, complexconc, &incx, &beta, indivconc, &incx); 
}
