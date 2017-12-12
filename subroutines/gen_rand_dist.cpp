#include "gen_rand_dist.h"
#include "rand_gsl.h"

void gen_rand_dist(int Nrand, double *cumul, double *abund, double *randvals, int Nconstrain, int *constrain, double *net_abund)
{
  int i;
  int s;
  double rnum;
  double df;
  double n1;
  int flag;
  int c;
  for(i=0;i<Nrand;i++){
    
    flag=0;
    for(c=0;c<Nconstrain;c++){
      if(i==constrain[c]){
	randvals[i]=net_abund[i];
	flag=1;
      }
    }
    if(flag==0){
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
  
}
