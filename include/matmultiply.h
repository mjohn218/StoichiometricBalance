#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#define dgemm dgemm_ 
#define dgemv dgemv_ 


extern "C" void dgemm_(char *transa, char *transb, int *m, int *
		       n, int *k, double *alpha, double *a, int *lda, 
		       double *b, int *ldb, double *beta, double *c, int 
		       *ldc);


extern "C" void dgemv_(char *trans, int *m, int *n, double *
		 alpha, double *a, int *lda, double *x, int *incx, 
		 double *beta, double *y, int *incy);
void matmultiply_nomkl(int nmol, int ncomplex, double *indivconc, double *complexconc, double *A);
/*this version below uses BLAS routines*/
void matmultiply(int nmol, int ncomplex, double *indivconc, double *complexconc, double *A);
