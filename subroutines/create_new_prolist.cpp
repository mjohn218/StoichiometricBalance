#include <vector>
#include <string>
#include "pro_classes.h"
#include "read_protein_class.h"
#include "construct_net.h"

void create_new_prolist( vector<ProteinClass> &proList, vector<ProteinClass> &proListKD,  vector<int> proKDID, int nProKD)
{
  ProteinClass temp(-1);
  int nConstrain=0;
  int t=0;
  int Npro=proList.size();
  for(int i=0;i<Npro;i++){
    bool keep=true;
    for(int p=0;p<nProKD;p++){
      if(proKDID[p]==i){
	//this protein is knocked down, exclude!
	keep=false;
	break;
      }
    }
    if(keep==true){
      temp.name=proList[i].name;
      temp.copies=proList[i].copies;
      temp.nInterfaces=0;//haven't added these yet, initialize to zero
      temp.isConstrained=proList[i].isConstrained;
      temp.index=t;
      /*Only populate these elements, will reconstruct interface network from scratch.*/
      proListKD.push_back(temp);//proList[i]);
      
      if(proList[i].isConstrained==true){
	nConstrain++;
	//constrainKD.push_back(t);
      }
      t++;//keeps track of number of proteins in proListKD
    }
  }
  
  cout <<" -------------------- "<<endl;
  cout <<" NUMBER OF PROTEINS IN KNOCKEDDOWN LIST: "<<proListKD.size()<<endl;
  cout <<" NUMBER OF CONSTRAINED PROTEIN ABUNDANCES in KnockDown : "<<nConstrain<<endl;



}
