#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_copt(int Nedge, double *complexopt, ofstream &outfile, int Nc, int *constrain)
{
  int i, j;
  outfile<<"Constraints: "<<'\t';
  for(i=0;i<Nc;i++)
    outfile<<constrain[i]<<'\t';
  outfile<<endl;

  for(i=0;i<Nedge;i++){
    outfile<<complexopt[i]<<endl;
    
  }
  
}
