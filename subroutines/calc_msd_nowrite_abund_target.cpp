#include "calc_abund_target_metrics.h"

double calc_msd_nowrite_abund_target(int Npro, double *proconc, double *abund)
{
  int i;
  double sum=0;
  double df, df2;

  for(i=0;i<Npro;i++){
    if(abund[i]>0){
      df=proconc[i]-abund[i];
      df2=df*df;
      sum+=df2;

    }    else{
      df=-1;

    }
    
  }
  double msd=sum*1.0/(1.0*Npro);
  return msd;
  
}
