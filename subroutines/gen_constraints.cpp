#include "gen_constraints.h"
#include "rand_gsl.h"

void gen_constraints(int Nconstrain, int *constrain, int Npro, double *abund)
{
  int i, j;
  double rnum;
  int ind;
  i=0;
  int sflag;
  while(i<Nconstrain){
    rnum=rand_gsl()*1.0*Npro;
    ind=int(rnum);
    while(abund[ind]<=0){
      rnum=rand_gsl()*1.0*Npro;
      ind=int(rnum);
    }
    constrain[i]=ind;
    if(i>0){
      //make sure you don't pick the same protein more than once
      sflag=0;
      for(j=0;j<i;j++){
	if(constrain[j]==constrain[i])
	  sflag=1;
      }
      if(sflag==1)
	i-=1;//reselect
    }
    i++;
  }

}
