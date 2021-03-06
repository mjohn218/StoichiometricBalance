#include "pro_classes.h"
#include "fitness_constrain_to_target.h"

double fitness_constrain_sqr2(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund)
{
  int i, j, k;
  int i1, i2, ni;
  double fit2=0;
  double chg;
  int p1;
  double c1;
  double dni;
  /*add in a constraint for the proteins that have a fixed concentration */
  for(i=0;i<Nconstrain;i++){
    p1=constrain[i];
    c1=net_abund[p1];
    dni=1.0*wholep[p1].ninterface;
    for(j=0;j<wholep[p1].ninterface;j++){
      i1=wholep[p1].valiface[j];
      chg=indivconc[i1]-c1;
      fit2+=chg*chg/dni;
    }
  }
  double npair;
  for(i=0;i<nwhole;i++){
    ni=wholep[i].ninterface;
    npair=ni*(ni-1)*1.0/2.0;
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      for(k=j+1;k<ni;k++){
	i2=wholep[i].valiface[k];
	chg=indivconc[i1]-indivconc[i2];
	//cout <<"prot: "<<i<<" if1: "<<i1<<" if2: "<<i2<<" chg: "<<chg<<endl;
	fit2+=chg*chg/(1.0*npair);
      }
    }
  }
  return fit2;
}  
