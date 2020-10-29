#include "build_Amatrix.h"


void build_Amatrix_with_Stoich(int nmol, double *A, int * numpartners, int **Speclist, double **Stoichlist)
{
  int i, j;
  int part;
  int t=0;
  double sval1, sval2;
  for(i=0;i<nmol;i++){
    for(j=0;j<numpartners[i];j++){
	part=Speclist[i][j];
	sval1=Stoichlist[i][j];
      
	if(part>i){
	//new reaction

	//rows of A are for each interface, columns are for each complex
	A[t*nmol+i]=sval1;
	for(int s=0;s<numpartners[part];s++){
	    if(Speclist[part][s]==i)
		A[t*nmol+part]=Stoichlist[part][s];//this should be overwritten by below.
	}
	t++;
      }
	
      if(part==i){
	A[t*nmol+i]=2;
	t++;
      }
    }
  }

  /* cout <<"-------------------"<<endl;
     cout <<" Write A MATRIX TRANSPOSE (First row is all interfaces in complex 1) to amatrixT.out file: "<<endl;
     ofstream amatrix("amatrixT.out");
     for(i=0;i<Nif*Ncomplex;i++){
     if((i+1)%Nif==0)
     amatrix <<A[i]<<endl;
     else
     amatrix <<A[i]<<'\t';
     }
     cout <<"-------------------"<<endl; 
  */
  
}
