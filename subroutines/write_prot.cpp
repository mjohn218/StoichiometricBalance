#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_prot(int nwhole, Protein *wholep, double *conc, ofstream &outfile, double &nprot )
{
  int i, j;
  int ni, iface;
  nprot=0;
  for(i=0;i<nwhole;i++){
    ni=wholep[i].ninterface;
    nprot+=conc[wholep[i].valiface[0]];
    for(j=0;j<ni;j++){
      iface=wholep[i].valiface[j];
      outfile<<conc[iface]<<'\t';
    }
    outfile<<endl;
  }
  
}
