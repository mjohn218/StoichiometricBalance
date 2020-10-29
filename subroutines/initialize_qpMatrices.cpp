#include <vector>
#include <string>
#include "pro_classes.h"
#include "constrainParms.h"
#include "initialize_qp.h"
/*
  Fill in all the matrices used in the QP calculation, based on the constrained protein copy numbers and the network structure, and stoichiometry.
*/
using namespace std;

void initialize_qpMatrices(int Ncomplex, int Nif, int nConstrain, constrainParms &plist, Protein *wholep, vector<int> constrain, double *A, int *p_home, double *xupp, char *ixupp, double *xlow, char *ixlow, int *irowQ, int *jcolQ,  double *clow, char *iclow, double *cupp, char *icupp, int *irowC, int *jcolC, double *dC, double *H, double *ZA)
{

  int nx=Ncomplex;
  int i, j;
  for(i=0;i<nx;i++){
    xupp[i]=0;
    ixupp[i]=0;
    xlow[i]=0;
    ixlow[i]=0;
  }
  int t=0;
  for(i=0;i<Ncomplex;i++){
    for(j=i;j<Ncomplex;j++){
      jcolQ[t]=i;
      irowQ[t]=j;
      //dQ[t]=Q[i*Ncomplex+j];
      t++;
    }
  }

  for(i=0;i<Ncomplex;i++){
    clow[i]=plist.min_complex;
    iclow[i]=1;
    cupp[i]=plist.max_complex;//In reality , this could be set for each element of the array, instead of using same value for all!
    if(plist.max_complex>0)
	icupp[i]=1;
    else
	icupp[i]=0;
  }
  
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }
  int p1,p2,  nint;
  double cof;
  for(i=0;i<Nif;i++){
    p1=p_home[i];
    nint=wholep[p1].ninterface;
    cof=1.0;
    H[i*Nif+i]=(nint-1)/cof;
    //cout <<"Hi: "<<H[i*Nif+i]<<" nint: "<<nint<<endl;
    
    for(j=i+1;j<Nif;j++){
      p2=p_home[j];
      if(p2==p1){
	H[i*Nif+j]=-1.0/cof;
	H[j*Nif+i]=-1.0/cof;
      }
      else{
	H[i*Nif+j]=0;
	H[j*Nif+i]=0;
	
      }
    }
  }
  double *Z=new double[Nif*Nif];
  //Start out Z as zero matrix
  for(i=0;i<Nif;i++){      
    Z[i*Nif+i]=0;
  }
  int ntmp;
  
  for(i=0;i<nConstrain;i++){
      for(j=0;j<wholep[constrain[i]].ninterface;j++){
	  ntmp=wholep[constrain[i]].valiface[j];
	  Z[ntmp*Nif+ntmp]=1;
      }
  }

  /*define Z*A matrix*/
  for(i=0;i<Nif;i++){
    if(Z[i*Nif+i]==1){
      //copy A into ZA
      for(j=0;j<Ncomplex;j++){
	ZA[j*Nif+i]=A[j*Nif+i];
      }
    }else{
      for(j=0;j<Ncomplex;j++)
	ZA[j*Nif+i]=0;
    }
  }

  delete[] Z;
}
