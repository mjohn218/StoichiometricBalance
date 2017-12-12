#include "pro_classes.h"
#include "constrainParms.h"
#include "read_proinput.h"

void read_protlist(int nwhole, Protein *wholep, int nmol, int *p_home, ifstream &protfile)
{
  //read in each proteins number of interfaces and their identities
  int i, j;
  int nint;
  int ig;
  int ct=0;
  int iface;
  for(i=0;i<nwhole;i++){
    protfile >>ig;//protein number, =i
    protfile >>nint;//number of protein i's interfaces
    wholep[i].ninterface=nint;
    cout <<"Protein: "<<i<<" Numinterfaces: "<<wholep[i].ninterface<<endl;
    for(j=0;j<nint;j++){
      protfile >>iface;
      wholep[i].valiface[j]=iface;
      p_home[iface]=i;
    }
    ct+=nint;
  }
  if(ct!=nmol){
    cerr<<"Number of interfaces assigned to proteins does not match network interface numbers!"<<endl;
    exit(1);
  }
  
}
