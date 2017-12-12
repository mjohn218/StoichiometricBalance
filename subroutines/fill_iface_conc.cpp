#include "pro_classes.h"
#include "fill_iface_conc.h"

void fill_iface_conc(int Npro, Protein *wholep, double *indivconc, double *proconc)
{
  int i, j, k;

  int i1, ni;
  
  for(i=0;i<Npro;i++){
    ni=wholep[i].ninterface;

    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      indivconc[i1]=proconc[i];
    }
  }
  
  
}
