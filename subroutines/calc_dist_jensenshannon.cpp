#include "calc_abund_target_metrics.h"

double calc_dist_jensenshannon(int Npro, double *proconc, double *abund, string *genid)
{
  /*See what is the discrepancy between predicted and actual expression levels, by
    turning both of them into probability distributions*/
  int i;
  double sum=0;
  double sumrev=0;
  double lg, revlg;
  double pred[Npro];
  double pobs[Npro];
  double N1=0;
  double N2=0;
  for(i=0;i<Npro;i++){
    if(abund[i]>0){
      N2+=abund[i];
      N1+=proconc[i];
    }
  }
  for(i=0;i<Npro;i++){
    if(abund[i]>0){
      pred[i]=proconc[i]*1.0/(N1*1.0);
      pobs[i]=abund[i]*1.0/(N2*1.0);
    }
  }
  //Now that they are converted to probability, need to get M=(P+Q)/2
  double M[Npro];
  for(i=0;i<Npro;i++)
    M[i]=0.5*(pred[i]+pobs[i]);
  
  //Get KL between P and M
  for(i=0;i<Npro;i++){
    if(abund[i]>0 && pred[i]>0){
      lg=pred[i]*log(pred[i]*1.0/(1.0*M[i]));
      //      revlg=pobs[i]*log(pobs[i]*1.0/(pred[i]*1.0));
      
      sum+=lg;
      //sumrev+=revlg;
    } 
    
  }
  //Get KL between Q and M
  for(i=0;i<Npro;i++){
    if(abund[i]>0 && pred[i]>0){
      //lg=pred[i]*log(pred[i]*1.0/(1.0*pobs[i]));
      lg=pobs[i]*log(pobs[i]*1.0/(M[i]*1.0));
      
      sum+=lg;
      //sumrev+=revlg;
    } 
    
  }
  //Divide by 2
  sum=0.5*sum;
  //Take sqrt to get distance
  double ent=sqrt(sum);
  return ent;
  
}
