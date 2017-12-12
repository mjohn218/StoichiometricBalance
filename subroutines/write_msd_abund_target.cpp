#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_msd_abund_target(int Npro, double *proconc, double *abund, ofstream &outfile, string *genid, double *randvals)
{
  /*See what is the discrepancy between predicted and actual expression levels*/
  int i;
  double sum=0;
  double df, df2;
  outfile<<"Pro"<<'\t'<<"Name" <<'\t'<<"ExpAbund"<<'\t'<<"Predict"<<'\t'<<"Random"<<'\t'<<"Diff"<<'\t'<<"Diff2"<<'\t'<<"SumDiff2"<<endl;
  for(i=0;i<Npro;i++){
    if(abund[i]>0){
      df=proconc[i]-abund[i];
      df2=df*df;
      sum+=df2;
      outfile<<i<<'\t'<<genid[i]<<'\t'<<abund[i]<<'\t'<<proconc[i]<<'\t'<<randvals[i]<<'\t'<<df<<'\t'<<df2<<'\t'<<sum<<endl;
    }    else{
      df=-1;
      outfile<<i<<'\t'<<genid[i]<<'\t'<<abund[i]<<'\t'<<proconc[i]<<'\t'<<randvals[i]<<'\t'<<"00"<<'\t'<<"00"<<'\t'<<sum<<endl;
    }
    
  }

  
}
