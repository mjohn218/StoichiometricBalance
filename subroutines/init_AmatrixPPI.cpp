#include "pro_classes.h"
#include "init_readref_Amatrix.h"

void init_AmatrixPPI(int Npro, int Nedge, int Nif, double *A, int *e1num, int *e2num, Protein *wholep, int *e1int, int *e2int)
{
  /*Gives each protein a new interface for each edge!*/
  int i, j;
  int part;
  int t=0;
  int myind[Npro];
  for(i=0;i<Npro;i++)
    myind[i]=0;
  
  int nc=0;
  int p1, p2;
  int i1, i2;

  for(i=0;i<Nedge;i++){
    p1=e1num[i];
    p2=e2num[i];
    i1=wholep[p1].valiface[myind[p1]];
    i2=wholep[p2].valiface[myind[p2]];
    e1int[i]=i1;
    e2int[i]=i2;
    //cout <<"protein 1: "<<p1<<" myind: "<<myind[p1]<<" mypartners: "<<wholep[p1].ninterface<<endl;
    //cout <<"protein 2: "<<p2<<" myind: "<<myind[p2]<<" mypartners: "<<wholep[p2].ninterface<<endl;
	
    if(p1==p2){
      i2=i1;
      A[i*Nif+i1]=2;
      myind[p1]++;
    }
    else{
      myind[p1]++;
      myind[p2]++;
      //in this case each interface has a single partner
      //rows are for each interface, columns are for each complex
      A[i*Nif+i1]=1;
      A[i*Nif+i2]=1;
    }
    nc++;
    //cout <<"ncomplex: "<<nc<<" p1: "<<p1<<" p2: "<<p2<<endl;
    
  }
   

}
