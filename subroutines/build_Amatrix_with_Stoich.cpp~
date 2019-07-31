#include "build_Amatrix.h"


void build_Amatrix(int nmol, double *A, int * numpartners, int **Speclist)
{
  int i, j;
  int part;
  int t=0;
  for(i=0;i<nmol;i++){
    for(j=0;j<numpartners[i];j++){
      part=Speclist[i][j];
      if(part>i){
	//new reaction

	//rows of A are for each interface, columns are for each complex
	A[t*nmol+i]=1;
	A[t*nmol+part]=1;
	t++;
      }
      if(part==i){
	A[t*nmol+i]=2;
	t++;
      }
    }
  }
}
