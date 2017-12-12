#include "mut_complexconc.h"
#include "utility_calls.h"

void mutate_glob(int ncomplex, int globswap, int minnum, double *complexmut)
{
  /*Mutate each complex a little*/
  int i;
  double rnum, dchg;
  int nchg;
  int tot=0;
  double xnew;
  for(i=0;i<ncomplex;i++){
    rnum=trand()*1.0/(1.0*RAND_MAX);
    /*choose a number of complexes to delete*/
    dchg=globswap*2.0*trand()/(RAND_MAX*1.0)-globswap;
    nchg=int(round(dchg));
    xnew=complexmut[i]+nchg;
    if(xnew<minnum)
      complexmut[i]=minnum;
    else
      complexmut[i]=complexmut[i]+nchg;

    
  }
  
}
