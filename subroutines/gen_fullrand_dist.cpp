#include "gen_fullrand_dist.h"
#include "rand_gsl.h"

void gen_fullrand_dist(int Nrand, double *cumul, double *abund, double *randvals)
{
  int i;
  int s;
  double rnum;
  double df;
  double n1;
  int flag;
  int c;
  for(i=0;i<Nrand;i++){

    rnum=rand_gsl()*1.0;
    /*reach rnum is a probability, will allow us to select from the
      distribution*/
    s=0;
    while(rnum>cumul[s]){
      s++;
    }
    /*now that we have the element, linearly interpolate*/
    df=(cumul[s]-rnum)/(cumul[s]-cumul[s-1]);
    n1=df*(abund[s]-abund[s-1])+abund[s-1];
    randvals[i]=n1;
  }
  
  
}
