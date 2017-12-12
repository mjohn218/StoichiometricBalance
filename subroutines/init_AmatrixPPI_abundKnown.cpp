#include "pro_classes.h"
#include "init_readref_Amatrix.h"

void init_AmatrixPPI_abundKnown(int Npro, int Nedge, int Nif, int *e1num, int *e2num, Protein *wholep, int *e1int, int *e2int, double *abund)
{
  /*formerly build_AmatrixPPI_both.
    this routine does not actually assign the Amatrix, and gives new interface
    to each edge only if the abundance of the protein is known.
  */
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
      if(abund[p1]<0)
	myind[p1]++;
      //      else{
	//this protein has a single interface,so don't advance it   
      //}
    }
    else{
      if(abund[p1]<0){
	myind[p1]++;
      }
      if(abund[p2]<0){
	myind[p2]++;
	//in this case each interface has a single partner
	//rows are for each interface, columns are for each complex
      }
    }
    nc++;
    //cout <<"ncomplex: "<<nc<<" p1: "<<p1<<" p2: "<<p2<<endl;
    
  }
   

}
