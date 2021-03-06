#include "pro_classes.h"
#include "functionCalls.h"


void get_variances(int Npro, int Nif, Protein *wholep, vector<double> R, double *ivars, double *imean)
{

  double sum1,sum2;
  int i,j;
  int i1,ni;

  //Calculate mean and variance
  
  
  for(i=0;i<Npro;i++){
    sum1=0;
    sum2=0;
    ni=wholep[i].ninterface;
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      //cout << i1 << ": " << R[i1] << endl;
      sum1+=R[i1];
      sum2+=R[i1]*R[i1];
    }
    //Now divide
    sum1/=(1.0*ni); //This is the mean
    sum2/=(1.0*ni); //This is E[x^2]
    ivars[i]=sum2 - sum1*sum1; //This is the variance
    imean[i]=sum1;
  }

}
