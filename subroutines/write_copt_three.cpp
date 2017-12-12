#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_copt_three(int Nedge, double *complexopt, double *restart, double *restart2, ofstream &outfile, int Nc, int *constrain )
{
  int i, j;
  outfile<<"Constraints: "<<'\t';
  for(i=0;i<Nc;i++)
    outfile<<constrain[i]<<'\t';
  outfile<<endl;
  for(i=0;i<Nedge;i++){
    outfile<<complexopt[i]<<'\t'<<restart[i]<<'\t'<<restart2[i]<<endl;
    
  }
  
}
