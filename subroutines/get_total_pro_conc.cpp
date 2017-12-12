#include "pro_classes.h"
#include "get_total_pro_conc.h"

double get_total_pro_conc(int nwhole, Protein *wholep, double *conc) 
{
  int i, j;
  int ni, iface;
  double nprot=0;
  for(i=0;i<nwhole;i++){
    nprot+=conc[wholep[i].valiface[0]];
  }
  return nprot;
}
