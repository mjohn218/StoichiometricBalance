#include "pro_classes.h"
#include "constrainParms.h"
#include "read_proinput.h"

void read_ref_ProIIN(int nwhole, Protein *wholep, int &nrefinterface, ifstream &reffile, int *npartners, int **Speclist, int *p_home)
{
  /*read in each proteins number of interfaces and their identities
    then read in IIN connectivities. 
  */
  int i, j;
  int nint;
  int ig;
  int ct=0;
  int iface;
  for(i=0;i<nwhole;i++){
    reffile >>ig;//protein number, =i
    reffile >>nint;//number of protein i's interfaces
    wholep[i].ninterface=nint;
    
    //cout <<"Protein: "<<i<<" Numinterfaces: "<<wholep[i].ninterface<<endl;
    if(nint==-1){
    }else{
      for(j=0;j<nint;j++){
	reffile >>iface;
	wholep[i].valiface[j]=iface;
	p_home[iface]=i;
      }
    }
    ct+=nint;
  }
  nrefinterface=ct;
  cout <<"N reference interfaces: "<<nrefinterface<<endl;
  for(i=0;i<nrefinterface;i++){
    reffile>>ig;
    reffile >>npartners[i];
    //cout <<ig<<'\t'<<npartners[i]<<'\t';
    for(j=0;j<npartners[i];j++){
      reffile>>Speclist[i][j];
      //cout <<Speclist[i][j]<<'\t';
    }
    //    cout <<endl;
  }
  
  reffile.close();

}
