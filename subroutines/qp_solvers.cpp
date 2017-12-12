#include "pro_classes.h"
#include "constrainParms.h"
#include "qp_solvers.h"
#include "matmultiply.h"
#include "gen_fullrand_dist.h"
#include "doubleLexSort.h"



void qp_init(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *randvals)
{

  int i, j;
  
  /*Initialize the Protein interface concentrations, from the yeast distribution*/
  gen_fullrand_dist(Nif, cumul, yeast_abund, randvals);
  //  cout <<"individual protein concentrations: "<<endl;
  numstart=0;
  for(i=0;i<Nif;i++)
    indivconc[i]=randvals[i];
  
  int i1, i2, ni;

  /*Use protein constraints as initial values for those protein interfaces concentrations*/
  int id;
  for(i=0;i<Nconstrain;i++){
    id=constrain[i];
    ni=wholep[id].ninterface;
    for(j=0;j<ni;j++){
      i1=wholep[id].valiface[j];
      indivconc[i1]=abund[id];
    }
  }
  for(i=0;i<Npro;i++){
    i1=wholep[i].valiface[0];
    numstart+=(int)indivconc[i1];
  }

  
  
  
  /*Now try to find an optimal value for the complex concentrations*/
  
  /*In QP, trying to minimize the objective function 1/2 x'Qx +c'x
    subject to constraints on x, either Cx>=d, (where C is Identity and
    D is zeros) or with bounds set as x>=0*/
  /*define the Q matrix, first just for optimizing Ax=b with x>0, by minimizing (Ax-b)^2*/
  /*Q matrix is 2*A'A, c' vector is -2*b'A, and RHS is -b'b constraint is for x>=0*/
  
  
  double Rhs=0.0;


  for(i=0;i<Nif;i++){
    //    cout <<"B[i]: "<<bvec[i]<<endl;
    Rhs+=indivconc[i]*indivconc[i];
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
  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &lda, indivconc, &incx, &beta, c,  &incx); 
  
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
  

  /*store lower triangle of Q in dQ, store matrix indices of those elements in i, j*/
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
  //  cout <<"Num nonzero elements in Q: "<<nnzQ<<" total elements: "<<Ncomplex*(Ncomplex+1)/2<<endl;
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
  
  const int mz   = Ncomplex;

  const int nnzC = Ncomplex;
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }
  

  QpGenSparseMa27 * qp 
    = new QpGenSparseMa27( nx, my, mz, nnzQ, nnzA, nnzC );
  

  QpGenData      * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
								      c,      irowQ,  nnzQ,   jcolQ,  dQ,
								      xlow,   ixlow,  xupp,   ixupp,
								      irowA,  nnzA,   jcolA,  dA,     b,
								      irowC,  nnzC,   jcolC,  dC,
								      clow,   iclow,  cupp,   icupp );
  

  QpGenVars      * vars 
    = (QpGenVars *) qp->makeVariables( prob );
  QpGenResiduals * resid 
    = (QpGenResiduals *) qp->makeResiduals( prob );
  

  GondzioSolver  * s     = new GondzioSolver( qp, prob );
  
  //  if( !quiet ) s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  
  if( ierr == 0 ) {
    //    cout.precision(4);
    cout << "Solution: \n";
    //    vars->x->writefToStream( cout, "x[%{index}] = %{value}" );
    vars->x->copyIntoArray(complexconc);
  } else {
    cout << "Could not solve the problem.\n";
    cout <<"generate random start point for complexes "<<endl;
    gen_fullrand_dist(Ncomplex, cumul, yeast_abund, complexconc);

    matmultiply(Nif,  Ncomplex, indivconc, complexconc, A);
    numstart=0;
    for(i=0;i<Npro;i++){
      i1=wholep[i].valiface[0];
      numstart+=(int)indivconc[i1];
    }
  }
  //      cout ierr;
  //double *xsol=new double[Ncomplex];
  double objective=prob->objectiveValue(vars)+ Rhs;
  /*Now solve Ax=b'*/
  trans='N';
  alpha=1.0;
  beta=0;
  incx=1;
 //  double *newsol=new double[Nif];
//   dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &Nif, complexconc, &incx, &beta, newsol, &incx); 
//   cout <<"new predicted bvector from optimal x vector: "<<endl;
//   ofstream bfile("Bpred.out");
//   bfile<<"Actual  Predicted"<<endl;
//   for(i=0;i<Nif;i++){
//     bfile<<indivconc[i]<<' '<<newsol[i]<<endl;
//   }
  
  delete qp;
  delete vars;
  delete prob;
  delete resid;
  delete s;
  delete dA;
  delete b;
  delete irowA;
  delete jcolA;
  
}
void qp_solve(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *H, double *ZA, double *Q2, double ascale, int *p_home)
{

  int i, j;

  cout << "Inside qp_solve" << endl;
  
  /*Now try to find an optimal value for the complex concentrations*/
  
  /*In QP, trying to minimize the objective function 1/2 x'Qx +c'x
    subject to constraints on x, either Cx>=d, (where C is Identity and
    D is zeros) or with bounds set as x>=0*/
  /*solve alpha(Ax-c)'Z(Ax-c)+(Ax)'H(Ax) where alpha is a constant scalar, c is the known concentrations, 
    Z picks out the known concentrations we're constraining, and H minimizes all interfaces on the same protein to
    zero, Z'=Z */
  /*Q matrix is 2*alpha*A'*Z*A+ 2*A'*H*A, c' vector is -2*alpha*c'*Z'*A, and RHS is -alpha*c'*Z*c constraint is for x>=0*/
  
  
  double Rhs=0.0;

  int p1;
  int ni;

  for(i=0;i<Nconstrain;i++){
    p1=constrain[i];
    ni=wholep[p1].ninterface;
    Rhs+=ascale*abund[p1]*abund[p1]*ni;//this is c'Z*c
  }
  int k, l;
   
  cout <<"RHS: "<<Rhs<<endl;
  
  /*A matrix is Ninf rows by Ncomplex columns*/
  char trans='T';
  //  double alpha=-2.0;
  int lda=Nif;
  int incx=1;
  double beta=0;
  
  /*define c vector*/
  //  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &lda, indivconc, &incx, &beta, c,  &incx); 
  for(j=0;j<Ncomplex;j++){
    c[j]=0;
    for(i=0;i<Nif;i++){
      //    do c'*Z*A
      p1=p_home[i];
      c[j]+=-2.0*ascale*abund[p1]*ZA[j*Nif+i];
      
    }
  }
 
  char transa='T';
  char transb='N';
  double alph2=2.0*ascale;
  double bet2=0.0;

  /*define Q matrix, 2 parts*/

  dgemm(&transa, &transb, &Ncomplex, &Ncomplex, &Nif, &alph2, A, &Nif, ZA, &Nif, &bet2, Q, &Ncomplex); 
  /*This one is A'*H*A*/
  
  double one=1.0;
  double two=2.0;
    
  
  dgemm(&transb, &transb, &Nif, &Ncomplex, &Nif, &one, H, &Nif, A, &Nif, &bet2, Q2, &Nif); 
  //A'*(HA)   HA is stored in Q2
  
  
  dgemm(&transa, &transb, &Ncomplex, &Ncomplex, &Nif, &two, A, &Nif, Q2, &Nif, &one, Q, &Ncomplex);

    
  /*Read Q, cvec, and zerovec into the QP program and
	it should produce the x-vector, or complexconc!*/
  int usage_ok = 1, quiet = 0;
  /*Data structures*/
  const int nx   = Ncomplex;
  //defined above:   double    c[]  = ;//{ 1.5,  -2 };
  

  /*store lower triangle of Q in dQ, store matrix indices of those elements in i, j*/
  //  const int nnzQ = Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2

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
  cout <<"Num nonzero elements in Q: "<<nnzQ<<" total elements: "<<Ncomplex*(Ncomplex+1)/2<<endl;
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
  
  const int mz   = Ncomplex;

  const int nnzC = Ncomplex;
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }


  QpGenSparseMa27 * qp 
    = new QpGenSparseMa27( nx, my, mz, nnzQ, nnzA, nnzC );
  

  QpGenData      * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
								      c,      irowQ,  nnzQ,   jcolQ,  dQ,
								      xlow,   ixlow,  xupp,   ixupp,
								      irowA,  nnzA,   jcolA,  dA,     b,
								      irowC,  nnzC,   jcolC,  dC,
								      clow,   iclow,  cupp,   icupp );
  

  QpGenVars      * vars 
    = (QpGenVars *) qp->makeVariables( prob );
  QpGenResiduals * resid 
    = (QpGenResiduals *) qp->makeResiduals( prob );
  

  GondzioSolver  * s     = new GondzioSolver( qp, prob );
    
  
  // if( !quiet ) s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  cout <<"Ierr: "<<ierr<<endl;
  if( ierr == 0 ) {
    //cout.precision(4);
    //cout << "Solution: \n";
    //vars->x->writefToStream( cout, "x[%{index}] = %{value}" );
    vars->x->copyIntoArray(complexconc);
  } else {
    cout << "Could not solve the problem.\n";
    
    cout <<"Take current optimum! "<<endl;
    //vars->x->writefToStream( cout, "x[%{index}] = %{value}" );
    vars->x->copyIntoArray(complexconc);
  
    //    cout <<"total number of proteins : "<<numstart<<endl;
  }
  //      cout ierr;
  //double *xsol=new double[Ncomplex];
  double ob=prob->objectiveValue(vars);
  double objective=ob+ Rhs;
  cout <<"objective function: "<<ob<<" Plus RHS: "<<objective <<endl;

  
  delete qp;
  delete vars;
  delete prob;
  delete resid;
  delete s;
  delete dA;
  delete b;
  delete irowA;
  delete jcolA;
  
}




void qp_solve1(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *randvals, double *H, double *ZA, double *Q2, double ascale, int *p_home, double avg_degree)
{

  int i, j;
  
  /*Now try to find an optimal value for the complex concentrations*/
  
  /*In QP, trying to minimize the objective function 1/2 x'Qx +c'x
    subject to constraints on x, either Cx>=d, (where C is Identity and
    D is zeros) or with bounds set as x>=0*/
  /*solve alpha(Ax-c)'Z(Ax-c)+(Ax)'H(Ax) where alpha is a constant scalar, c is the known concentrations, 
    Z picks out the known concentrations we're constraining, and H minimizes all interfaces on the same protein to
    zero*/
  /*Q matrix is 2*alpha*A'*Z*A+ 2*A'*H*A, c' vector is -2*alpha*c'*Z'*A, and RHS is -alpha*c'*Z*c' constraint is for x>=0*/
  
  
  double Rhs=0.0;

  int p1;
  int ni;
  for(i=0;i<Nconstrain;i++){
    p1=constrain[i];
    ni=wholep[p1].ninterface;
    Rhs+=ascale*abund[p1]*abund[p1]*ni;//this is c'Z*c
  }
  int k, l;
  
  //  cout <<"RHS: "<<Rhs<<endl;
  
  /*A matrix is Ninf rows by Ncomplex columns*/
  char trans='T';
  //  double alpha=-2.0;
  int lda=Nif;
  int incx=1;
  double beta=0;
  
  /*define c vector*/
  //  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &lda, indivconc, &incx, &beta, c,  &incx); 
  for(j=0;j<Ncomplex;j++){
    c[j]=0;
    for(i=0;i<Nif;i++){
      //    do c'*Z*A
      p1=p_home[i];
      c[j]+=-2.0*ascale*abund[p1]*ZA[j*Nif+i];
      
    }
  }
  char transa='T';
  char transb='N';
  double alph2=2.0*ascale;
  double bet2=0.0;

  /*define Q matrix, 2 parts*/

  dgemm(&transa, &transb, &Ncomplex, &Ncomplex, &Nif, &alph2, A, &Nif, ZA, &Nif, &bet2, Q, &Ncomplex); 
  /*This one is A'*H*A*/
  double one=1.0;
  double two=2.0;
  
  dgemm(&transb, &transb, &Nif, &Ncomplex, &Nif, &one, H, &Nif, A, &Nif, &bet2, Q2, &Nif); 
  //A'*(HA)   HA is stored in Q2
  
  dgemm(&transa, &transb, &Ncomplex, &Ncomplex, &Nif, &two, A, &Nif, Q2, &Nif, &one, Q, &Ncomplex); 


  
  /*Read Q, cvec, and zerovec into the QP program and
	it should produce the x-vector, or complexconc!*/
  int usage_ok = 1, quiet = 0;
  /*Data structures*/
  const int nx   = Ncomplex;
  //defined above:   double    c[]  = ;//{ 1.5,  -2 };
  
  /*THIS VERSION DIFFERS FROM ABOVE IN DEFINTION JUST BELOW OF DQ[T], ABOVE VERSION CHECKS IF Q=0.*/
  /*store lower triangle of Q in dQ, store matrix indices of those elements in i, j*/
  const int nnzQ = Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2

  int t=0;
  for(i=0;i<Ncomplex;i++){
    for(j=i;j<Ncomplex;j++){
      jcolQ[t]=i;
      irowQ[t]=j;
      dQ[t]=Q[i*Ncomplex+j];
      t++;
    }
  }
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
  
  const int mz   = Ncomplex;

  const int nnzC = Ncomplex;
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }
  

  QpGenSparseMa27 * qp 
    = new QpGenSparseMa27( nx, my, mz, nnzQ, nnzA, nnzC );
  

  QpGenData      * prob = (QpGenData * ) qp->copyDataFromSparseTriple(
								      c,      irowQ,  nnzQ,   jcolQ,  dQ,
								      xlow,   ixlow,  xupp,   ixupp,
								      irowA,  nnzA,   jcolA,  dA,     b,
								      irowC,  nnzC,   jcolC,  dC,
								      clow,   iclow,  cupp,   icupp );
  

  QpGenVars      * vars 
    = (QpGenVars *) qp->makeVariables( prob );
  QpGenResiduals * resid 
    = (QpGenResiduals *) qp->makeResiduals( prob );
  

  GondzioSolver  * s     = new GondzioSolver( qp, prob );
  
  //  if( !quiet ) s->monitorSelf();
  int ierr = s->solve(prob,vars, resid);
  
  if( ierr == 0 ) {
    //    cout.precision(4);
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
    numstart=0;
    int i1;
    for(i=0;i<Npro;i++){
      i1=wholep[i].valiface[0];
      numstart+=(int)indivconc[i1];
    }
  }
  //      cout ierr;
  //double *xsol=new double[Ncomplex];
  double ob=prob->objectiveValue(vars);
  double objective=ob+ Rhs;
  cout <<"objective function: "<<ob<<" Plus RHS: "<<objective <<endl;
  /*Now solve Ax=b'*/
//   trans='N';
//   double alpha=1.0;
//   beta=0;
//   incx=1;
//   double *newsol=new double[Nif];
//  dgemv(&trans, &Nif, &Ncomplex, &alpha, A, &Nif, complexconc, &incx, &beta, newsol, &incx); 
  //  cout <<"new predicted bvector from optimal x vector: "<<endl;
//   ofstream bfile("Bpred.out");
//   bfile<<"Actual  Predicted"<<endl;
//   for(i=0;i<Nif;i++){
//     bfile<<newsol[i]<<endl;
//   }
  
  delete qp;
  delete vars;
  delete prob;
  delete resid;
  delete s;
  delete dA;
  delete b;
  delete irowA;
  delete jcolA;
  
}
