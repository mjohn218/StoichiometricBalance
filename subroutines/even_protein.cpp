#include "pro_classes.h"
#include "even_protein.h"

void even_protein(int nwhole, Protein *wholep, double *indivconc, double &nprot)
{
  int i, j, k;

  int i1, ni;
  double avg=0;
  nprot=0;
  double nav;
  for(i=0;i<nwhole;i++){
    ni=wholep[i].ninterface;
    avg=0;
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      avg+=indivconc[i1];
    }
    nav=avg/(ni*1.0);
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      indivconc[i1]=nav;
    }
    nprot+=nav;
    
  }
  
  
}
