#include <vector>
#include <string>
#include "pro_classes.h"
#include "read_protein_class.h"
#include "read_proinput.h"

void read_knockdowns(ifstream &kdfile, vector<string> &knockDowns, int Npro, vector<ProteinClass> proList, vector<int> &proKDID, int &nKD)
{
  
  cout<<"----------"<<endl;
  
  string pro1;
  ProteinClass temp(-1);
  while(!kdfile.eof()){
    //read in the names of each protein to knock down
    kdfile >>pro1;
    knockDowns.push_back(pro1);
  }
  //remove any duplicates
  nKD=knockDowns.size();
  for(int i=0;i<nKD;i++){
    for(int j=i+1;j<nKD;j++){
      if(knockDowns[i]==knockDowns[j]){
	cout <<"found duplicate KD, remove "<<endl;
	knockDowns[j]=knockDowns[nKD-1];//replace with last element
	nKD--;
      }
    }
  }
  /*Now compare each protein in the KD list with the full list of proteins.*/
  for(int i=0;i<nKD;i++){
    bool isFound=false;
    for(int j=0;j<Npro;j++){
      if(knockDowns[i]==proList[j].name){
	isFound=true;
	proKDID.push_back(j);
	cout <<" KNOCKDOWN PROTEIN: "<<knockDowns[i]<<endl;
	break;
      }
    }
    if(isFound==false){
      cout <<"CANNOT KNOCKDOWN protein: "<<knockDowns[i]<<" Not found in protein list, SKIPPING... "<<endl;
      knockDowns[i]=knockDowns[nKD-1];
      nKD--;
    }
  }



}
