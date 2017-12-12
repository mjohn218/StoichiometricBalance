#include "mut_complexconc.h"
#include "utility_calls.h"

void mutate_single(int ncomplex, int maxswap, int minnum, double *complexmut)
{
  /*mutate the number of protein complexes,
   either swap back and forth complexes (even numbers!)
   or mutate each complex a little*/
  double rnum=rand()*1.0/(1.0*RAND_MAX);
  /*choose a complex to delete from, and another to add to*/
  int c1=int(rnum*ncomplex);
  /*choose a number of complexes to delete*/
  double dchg=maxswap*2.0*rand()/(RAND_MAX*1.0)-maxswap;
  int nchg=int(round(dchg));
  
  double xnew=complexmut[c1]+nchg;//nchg can be positive or negative
  
  if(xnew>minnum){
    complexmut[c1]=complexmut[c1]+nchg;
  }
  //otherwise we skip this move
}
