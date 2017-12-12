#include "mut_complexconc.h"
#include "utility_calls.h"

void mutate_comp(int ncomplex, int maxswap, int minnum, double *complexmut)
{
  /*mutate the number of protein complexes,
   either swap back and forth complexes (even numbers!)
   or mutate each complex a little*/
  double rnum=rand()*1.0/(1.0*RAND_MAX);
  /*choose a complex to delete from, and another to add to*/
  int c1=int(rnum*ncomplex);
  int c2=int(1.0*rand()*ncomplex/(RAND_MAX*1.0));
  while(c2==c1){
    c2=int(1.0*rand()*ncomplex/(RAND_MAX*1.0));
  }
  /*choose a number of complexes to delete*/
  double dchg=maxswap*1.0*rand()/(RAND_MAX*1.0);
  int nchg=int(round(dchg));
  
  double xnew=complexmut[c1]-nchg;
  double x2=complexmut[c2]-nchg;
  if(xnew<minnum){
    //try subtracting off of c2 instead, otherwise skip the move
    if(x2>minnum){
      complexmut[c2]=complexmut[c2]-nchg;
      complexmut[c1]=complexmut[c1]+nchg;
    }
  }else{
    complexmut[c1]=complexmut[c1]-nchg;
    complexmut[c2]=complexmut[c2]+nchg;
  }

}
