#include "pro_classes.h"
#include "median_protein.h"
#include "utility_calls.h"

void median_protein(int nwhole, Protein *wholep, double *indivconc, double *proconc)
{
  int i, j, k;

  int i1, ni;
  double avg=0;
  long unsigned int index[MAXIFACE+1];
  double value[MAXIFACE+1];

  int id1, id2;
  double median;
  int mid;
  
  for(i=0;i<nwhole;i++){
    ni=wholep[i].ninterface;
    avg=0;
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      value[j+1]=indivconc[i1];
    }
    indexx(ni, &value[0], &index[0]);
    //cout <<"Protein: "<<i<<" ninterface: "<<ni<<endl;
    // for(j=0;j<ni;j++){
//       id1=index[j+1];
//       cout <<value[id1]<<endl;
//     }
    if(ni%2==0){
      mid=ni/2;
      id1=index[mid];
      id2=index[mid+1];
      median=0.5*value[id1]+0.5*value[id2];
    }else{
      mid=int(ni/2)+1;
      id1=index[mid];
      median=value[id1];
    }
    //cout <<"id1: "<<id1<<' '<<" mid: "<<mid<<endl;
    proconc[i]=median;
  }
  
  
}
