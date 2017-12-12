#include "pro_classes.h"
#include "qp_initconc.h"
#include "doubleLexSort.h"
#include "matmultiply.h"
#include "gen_fullrand_dist.h"


void qp_initconc(int Npro, int Nif, int Ncomplex, Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, double avg_degree)
{

  int i, j;
  
  /*Now try to find an optimal value for the complex concentrations. Protein interface concentrations are same for each protein,
   */
  
  /*In QP, trying to minimize the objective function 1/2 x'Qx +c'x
    subject to constraints on x, either Cx>=d, (where C is Identity and
    D is zeros) or with bounds set as x>=0*/
  /*define the Q matrix, first just for optimizing Ax=b with x>0, by minimizing (Ax-b)^2*/
  /*Q matrix is 2*A'A, c' vector is -2*b'A, and RHS is -b'b constraint is for x>=0*/
  
  double *Q=new double[Ncomplex*Ncomplex];
  double *c=new double[Ncomplex];
  double *bvec=new double[Nif];
  double Rhs=0.0;


  for(i=0;i<Nif;i++){
    bvec[i]=indivconc[i];
    //    cout <<"B[i]: "<<bvec[i]<<endl;
    Rhs+=-bvec[i]*bvec[i];
  }
  int k, l;
  
  //  cout <<"RHS: "<<Rhs<<endl;
  
  /*A matrix is Ninf rows by Ncomplex columns*/
  char trans='T';
  double alpha=-2.0;
  int lda=Nif;
  int incx=1;
  double beta=0;
  
  /*define c vector*/
  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &lda, bvec, &incx, &beta, c,  &incx); 
  
  char transa='T';
  char transb='N';
  double alph2=2.0;
  double bet2=0.0;
  /*define Q matrix*/
  dgemm(&transa, &transb, &Ncomplex, &Ncomplex, &Nif, &alph2, A, &Nif, A, &Nif, &bet2, Q, &Ncomplex); 
  
  /*Read Q, cvec, and zerovec into the QP program and
	it should produce the x-vector, or complexconc!*/
  int usage_ok = 1, quiet = 0;
  /*Data structures*/
  const int nx   = Ncomplex;
  //defined above:   double    c[]  = ;//{ 1.5,  -2 };
  
  /*bounds, should be zero*/
  double  *xupp=new double[nx];//0;//[] = { 20,   0 };
  char   *ixupp=new char[nx];//0;//[] = {  1,   0 };
  
  double  *xlow=new double[nx];//0;//[] = {  0,   0 };
  char   *ixlow=new char[nx];//0;//[] = {  1,   1 };
  for(i=0;i<nx;i++){
    xupp[i]=0;
    ixupp[i]=0;
    xlow[i]=0;
    ixlow[i]=0;
  }
  /*store lower triangle of Q in dQ, store matrix indices of those elements in i, j*/
  int tmpsize=Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2

  
  int    *irowQ=new int[tmpsize];// = {  0,   1,   1 }; 
  int    *jcolQ=new int[tmpsize];//[] = {  0,   0,   1 };
  double   * dQ=new double[tmpsize];//[] = {  8,   2,  10 };
  int t=0;
  
  for(i=0;i<Ncomplex;i++){
    for(j=i;j<Ncomplex;j++){
      if(Q[i*Ncomplex+j]!=0){
	jcolQ[t]=i;
	irowQ[t]=j;
	dQ[t]=Q[i*Ncomplex+j];
	t++;
      }
    }
  }
  //  const int nnzQ = Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2
  const int nnzQ=t;

  doubleLexSort(irowQ, nnzQ, jcolQ, dQ);
  //  cout <<"total number of lower triangle elements: "<<t<<" calculated: "<<nnzQ<<endl;
  //for(i=0;i<nnzQ;i++){
  // cout <<irowQ[i]<<' '<<jcolQ[i]<<' '<<dQ[i]<<endl;
  //}
  /*Equality matrix, nulled */
  int my         = 0;
  double * b     = 0;
  int nnzA       = 0;
  int * irowA    = 0;
  int * jcolA    = 0;
  double * dA    = 0;
  
  /*Inequality constraint*/
  /*      int mz=0;
	  double *clow=0;
	  double *cupp=0;
	  char *iclow=0;
	  char *icupp=0;
	  int nnzC=0;
	  int *irowC=0;
	  int *jcolC=0;
	  double *dC=0;*/
  
  const int mz   = Ncomplex;
  double *clow=new double[mz];
  char  *iclow=new char[mz];
  
  double *cupp=new double[mz];
  char  *icupp=new char[mz];
  for(i=0;i<mz;i++){
    clow[i]=0;
    iclow[i]=1;
    cupp[i]=0;
    icupp[i]=0;
  }
  const int nnzC = Ncomplex;
  int   *irowC=new int[nnzC];
  int   *jcolC=new int[nnzC];
  double   *dC=new double[nnzC];
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }
  
  cout <<"done allocating " <<endl; 
  QpGenSparseMa27 * qp 
    = new QpGenSparseMa27( nx, my, mz, nnzQ, nnzA, nnzC );
  
  cout <<"done sparse Ma27 "<<endl;
  QpGenData      * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
								      c,      irowQ,  nnzQ,   jcolQ,  dQ,
								      xlow,   ixlow,  xupp,   ixupp,
								      irowA,  nnzA,   jcolA,  dA,     b,
								      irowC,  nnzC,   jcolC,  dC,
								      clow,   iclow,  cupp,   icupp );
  
  cout <<"defined prob data structure "<<endl;
  QpGenVars      * vars 
    = (QpGenVars *) qp->makeVariables( prob );
  QpGenResiduals * resid 
    = (QpGenResiduals *) qp->makeResiduals( prob );
  
  cout<<" now solve: "<<endl;
  GondzioSolver  * s     = new GondzioSolver( qp, prob );
  
  if( !quiet ) s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  
  if( ierr == 0 ) {
    cout.precision(4);
    //    cout << "Solution: \n";
    //vars->x->writefToStream( cout, "x[%{index}] = %{value}" );
    vars->x->copyIntoArray(complexconc);
  } else {
    cout << "Could not solve the problem.\n";
    cout <<"generate random start point for complexes "<<endl;
    gen_fullrand_dist(Ncomplex, cumul, yeast_abund, complexconc);
    for(i=0;i<Ncomplex;i++)
      complexconc[i]/=avg_degree;
    matmultiply(Nif,  Ncomplex, indivconc, complexconc, A);
    
  }
  //      cout ierr;
  //double *xsol=new double[Ncomplex];

  /*Now solve Ax=b'*/
  trans='N';
  alpha=1.0;
  beta=0;
  incx=1;
  double *newsol=new double[Nif];
  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &Nif, complexconc, &incx, &beta, newsol, &incx); 
  cout <<"new predicted bvector from optimal x vector: "<<endl;
  //   ofstream bfile("Bpred.out");
  cout<<"Initial,   QP Predicted"<<endl;
  for(i=0;i<Nif;i++){
    cout<<indivconc[i]<<' '<<newsol[i]<<endl;
  }
  
  
}

