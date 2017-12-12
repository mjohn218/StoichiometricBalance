#include "accept_conc_mut.h"

void accept_conc_mut(int ncomplex, double &oldfit, double newfit, double &globopt, double *complexconc, double *complexmut, double *complexopt)
{
  int i;
  for(i=0;i<ncomplex;i++)
    complexconc[i]=complexmut[i];
  oldfit=newfit;
  if(newfit<globopt){
    globopt=newfit;
    for(i=0;i<ncomplex;i++)
      complexopt[i]=complexmut[i];
  }
}
