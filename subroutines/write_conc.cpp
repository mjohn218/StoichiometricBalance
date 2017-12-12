#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_conc(int nmol, double *conc, ofstream &outfile)
{
  int i;
  for(i=0;i<nmol;i++)
    outfile<<conc[i]<<endl;
  
}
