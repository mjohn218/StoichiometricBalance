#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_complex(int Nedge, double *conc, int *e1num, int *e2num, int *e1int, int *e2int, ofstream &outfile)
{
  int i;
  int p1, p2;
  outfile<<"pro1"<<'\t'<<"iface1"<<'\t'<<"pro2"<<'\t'<<"iface2"<<'\t'<<"abund"<<endl;
  for(i=0;i<Nedge;i++){
    p1=e1num[i];
    p2=e2num[i];
    outfile<<p1<<'\t'<< e1int[i]<<'\t'<<p2<<'\t'<<e2int[i]<<'\t'<<conc[i]<<endl;
  }
}
