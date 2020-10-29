#include <vector>
#include "ProteinClass.h"
#include "pro_classes.h"

#include "constrainParms.h"
#include "construct_net.h"


void init_matrix_file(ofstream &matrixFile, ofstream &matrixFileRB, vector<ProteinClass> proList, vector<double> proBalanced)
{
  

  
  int p1, i1;
  int Nif=0;
  double rb;
  
  int Npro=proList.size();
  vector<double> balancedPerPro;//populate this for all proteins used in SB.
  int n=0;
  vector<double> rbRatio;//only populate rbRatio if real copies>0.
  matrixFile<<" -- ";  
  for(p1=0;p1<proList.size();p1++){
    /*if p1 is in the KD list, need to skip it.*/
    //if(proList[p1].InterfaceList.size()>1){
    matrixFile<<'\t'<<proList[p1].name;
  }
  matrixFile<<endl;
  matrixFile<<"Orig";
  int index;
  for(p1=0;p1<proList.size();p1++){
    
    index=proList[p1].InterfaceList[0].globalIndex;
    balancedPerPro.push_back(proBalanced[index]);
    matrixFile<<'\t'<<proBalanced[index];
  }
  matrixFile<<endl;
  matrixFileRB<<" --";
  for(p1=0;p1<proList.size();p1++){
    if(proList[p1].copies>0){
      matrixFileRB<<'\t'<<proList[p1].name;
    }
  }
  matrixFileRB<<endl;
  matrixFileRB<<"Orig";
  for(p1=0;p1<proList.size();p1++){
    if(proList[p1].copies>0){
      rb=proList[p1].copies/balancedPerPro[p1];
      rbRatio.push_back(rb);
      matrixFileRB<<'\t'<<rb;
    }
  } 
  matrixFileRB<<endl;


  
}
