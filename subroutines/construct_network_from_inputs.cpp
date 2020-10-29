#include <vector>
#include "ProteinClass.h"
#include "pro_classes.h"

#include "constrainParms.h"
#include "construct_net.h"


/*Based on the proteins read in , and the list of edges read in, populate the class proList, where each edge is matched to proteins read in via the input files.
  This will thus ignore any edges read in that include proteins not listed in our list.
  
*/
void construct_network_from_inputs(int NedgeRead, int &nTotalInterfaces, int &maxInterfacePartners, int &nTotalEdges, vector<ProteinClass> &proList, vector<string> &edge1Interface, vector<string> &edge2Interface, vector<string> &edge1Protein, vector<string> &edge2Protein, vector<string> &stoichiometry, vector<string> filenames)
{

  Interface itemp(0);
  cout <<" -------------------- "<<endl;
  cout <<" MATCH EDGES BASED ON PROTEIN LIST FROM FILE "<<filenames[1]<<endl;

  double *stoichvalue1=new double[NedgeRead];
  double *stoichvalue2=new double[NedgeRead];
  string delimiter=":";
  int p;
  int p1, p2;
  int j;
  for(int i=0;i<NedgeRead;i++){
    /*its possible there is an edge in the network with a node not included */
    p1=-1;
    p2=-1;
    for(j=0;j<proList.size();j++){
      if(edge1Protein[i]==proList[j].name)
	p1=j;
      if(edge2Protein[i]==proList[j].name)
	p2=j;
    }
    if(p1==-1 || p2==-1){
	cout <<"CANNOT FIND BOTH PROTEINS IN EDGE: "<<i<<" skipping... " <<edge1Protein[i]<<' '<<edge2Protein[i]<<endl;
    }else{
	nTotalEdges++;
	//cout <<"Current edge: "<<i<<" proteins: "<<p1<<' '<<p2<<' '<<endl;
	bool doesExist1=false;//does the interface exist yet?
	bool doesExist2=false;//does the interface exist yet?
	/*Use this protein interaction to construct the interface interaction network*/
	int currIndex;
	int localIndex;
	//cout <<"P1: "<<p1<<" Current number of interfaces: "<<proList[p1].InterfaceList.size()<<" Interface to match: "<<edge1Interface[i]<<" edge iter: "<<nTotalEdges-1<<endl;
	for(int m=0;m<proList[p1].InterfaceList.size();m++){
	    //does this interface exist? Otherwise add it.
	    if(proList[p1].InterfaceList[m].name==edge1Interface[i]){
		//this interface already exists
		doesExist1=true;
		//cout <<"This interface already exists! "<<edge1Interface[i]<<" at edge: "<<nTotalEdges-1<<endl;
		localIndex=m;
		currIndex=proList[p1].InterfaceList[m].globalIndex;//which interface this is of all.
		
	    }
	}//looping over existing nInterfaces
	/*split the stoichiometry ratio into two values, stoichvalue1 and stoichvalue2*/
	
	string str=stoichiometry[i].substr(0, stoichiometry[i].find(delimiter));
	
	int start=stoichiometry[i].find(delimiter)+delimiter.length();
	stoichvalue1[i]=stod(str);
	string str2=stoichiometry[i].substr(start, stoichiometry[i].length());
	stoichvalue2[i]=stod(str2);
	//cout <<" first ratio value: "<<str<<" second value "<<str2<<" doubles: "<<stoichvalue1[i]<<' '<<stoichvalue2[i]<<endl;
	if(doesExist1==false){
	    
	    proList[p1].nInterfaces++;
	    localIndex=proList[p1].InterfaceList.size();//new interface on the protein, size before increasing
	    //create a new interface for this protein.
	    proList[p1].InterfaceList.push_back(itemp);//add this interface to the list for this protein.
	    proList[p1].InterfaceList[localIndex].name=edge1Interface[i];
	    //proList[p1].InterfaceList[localIndex].proteinHome=edge1Protein[i];
	    proList[p1].InterfaceList[localIndex].globalIndex=nTotalInterfaces;//index in full interface list
	    currIndex=nTotalInterfaces;
	   
	    nTotalInterfaces++;
	    //proList[p1].iFaceNames.push_back(edge1Interface[i]);//add this interface to the protein.
	    
	    
	}
	/**/
	proList[p1].InterfaceList[localIndex].proPartnerName.push_back(edge2Protein[i]);//add protein
	proList[p1].InterfaceList[localIndex].proPartnerIndex.push_back(p2);//add protein
	proList[p1].InterfaceList[localIndex].interfacePartnerName.push_back(edge2Interface[i]);//add interface
	proList[p1].InterfaceList[localIndex].stoichiometry.push_back(stoichvalue1[i]);
	int currsize=proList[p1].InterfaceList[localIndex].proPartnerName.size();
	if(currsize>maxInterfacePartners){
	    maxInterfacePartners=currsize;
	}

	/*Now do P2, add this interface and the interaction to the proteins*/
	//	cout <<"P2: "<<p2<<" Current number of interfaces: "<<proList[p2].InterfaceList.size()<<" Interface to match: "<<edge2Interface[i]<<endl;
	if(p1==p2 && edge1Interface[i]==edge2Interface[i]){
	    //this is a self interaction, don't add it again. 
	    //cout <<"SELF BINDING INTERFACE: "<<edge1Interface[i]<<' '<<edge2Interface[i]<<" on Protein: "<<p1<<endl;
	}else{
	    for(int m=0;m<proList[p2].InterfaceList.size();m++){
		//does this interface exist? Otherwise add it.
		if(proList[p2].InterfaceList[m].name==edge2Interface[i]){
		    //this interface already exists
		    doesExist2=true;
		    localIndex=m;
		    
		    currIndex=proList[p2].InterfaceList[m].globalIndex;//which interface this is of all.
		}
	    }//looping over existing nInterfaces
	    
	    
	    if(doesExist2==false){
		
		proList[p2].nInterfaces++;
		localIndex=proList[p2].InterfaceList.size();//new interface on the protein, current index
		//create a new interface for this protein.
		proList[p2].InterfaceList.push_back(itemp);//add this interface to the list for this protein.
		
		proList[p2].InterfaceList[localIndex].name=edge2Interface[i];
		//proList[p2].InterfaceList[localIndex].proteinHome=edge1Protein[i];
		proList[p2].InterfaceList[localIndex].globalIndex=nTotalInterfaces;//index in full interface list
		currIndex=nTotalInterfaces;
		nTotalInterfaces++;
		
	    }
	    /**/
	    proList[p2].InterfaceList[localIndex].proPartnerName.push_back(edge1Protein[i]);//add protein
	    proList[p2].InterfaceList[localIndex].proPartnerIndex.push_back(p1);//add protein
	    proList[p2].InterfaceList[localIndex].interfacePartnerName.push_back(edge1Interface[i]);//add interface
	    proList[p2].InterfaceList[localIndex].stoichiometry.push_back(stoichvalue2[i]);//add interface
	    currsize=proList[p2].InterfaceList[localIndex].proPartnerName.size();
	    
	    if(currsize>maxInterfacePartners){
		maxInterfacePartners=currsize;
	    }
	}//non self, add P2 interface
	
    }
  }//done looping over all edges in the input file.
  
  int i1;
  int nint;
  /*Now figure out globalIndex of each interface's partner interface*/
  for(p1=0;p1<proList.size();p1++){
    for(j=0;j<proList[p1].InterfaceList.size();j++){
      i1=proList[p1].InterfaceList[j].globalIndex;
      for(int k=0;k<proList[p1].InterfaceList[j].interfacePartnerName.size();k++){
	string i2name=proList[p1].InterfaceList[j].interfacePartnerName[k];
	int p2=proList[p1].InterfaceList[j].proPartnerIndex[k];
	//find the global index of this interface i2name on protein p2.
	int i2=-1;
	for(int s=0;s<proList[p2].InterfaceList.size();s++){
	  if(proList[p2].InterfaceList[s].name==i2name)
	    i2=proList[p2].InterfaceList[s].globalIndex;
	}
	proList[p1].InterfaceList[j].interfacePartnerIndex.push_back(i2);//add this to the interfaceList of partners
      }
    }
  }//found index of each interface's partner interface.

  /*NOW REMOVE ANY PROTEINS THAT WERE READ IN, BUT ARE NOT IN THE NETWORK*/
  int size=proList.size();
  for(int i=0;i<proList.size();i++){
    if(proList[i].InterfaceList.size()==0){
      cout <<"PROTEIN : "<<proList[i].name<<" HAS NO EDGES, REMOVE FROM LIST "<<endl;
      proList.erase(proList.begin()+i);
      i--;//erased a protein, so go back in case the next protein needs erasing. 
    }
  }
  //AS A CHECK, PRINT ALL PROTEINS IN ORDER:
  for(int i=0;i<proList.size();i++){
    cout<<"index: "<<i<<" name: "<<proList[i].name<<" isConstrain: "<<proList[i].isConstrained<<endl;
  }
  
  cout <<" ------------------- "<<endl;
  cout <<" TOTAL PROTEINS: "<<proList.size()<<" TOTAL INTERFACES: "<<nTotalInterfaces<<" MAX PARNTER PER INTERFACE: "<<maxInterfacePartners<<endl;
  cout <<" ---------------------" <<endl;
  cout <<" Total IIN Edges for these proteins: "<<nTotalEdges<<". Skipping a total of: "<<NedgeRead-nTotalEdges<<" edges from the file "<<filenames[2]<<" that did not have proteins listed in the file "<<filenames[1]<<endl;

  
  cout <<"----------------------------------"<<endl;
  cout <<" PRINT OUT PROTEIN INTERFACES and PROPERTIES : "<<endl;
  cout <<"----------------------------------"<<endl;
  /*Print out interaction network and Protein properties*/
  for(p1=0;p1<proList.size();p1++){
      cout <<" Protein index: "<<p1<<" Name: "<<proList[p1].name<<" copies: "<<proList[p1].copies<<" Total interfaces: "<<proList[p1].InterfaceList.size()<<" Name, index: "<<'\t';
      for(j=0;j<proList[p1].InterfaceList.size();j++){
	  i1=proList[p1].InterfaceList[j].globalIndex;
	  cout <<proList[p1].InterfaceList[j].name<<","<<proList[p1].InterfaceList[j].globalIndex<<'\t';
      }
      cout <<endl;
      /*Below commented out, printing all interface binding partners.*/
      /*
      for(j=0;j<proList[p1].InterfaceList.size();j++){
	  cout <<proList[p1].name<<'\t'<<proList[p1].InterfaceList[j].name<<" Partners: "<<proList[p1].InterfaceList[j].interfacePartnerName.size()<<endl;
	  for(int k=0;k<proList[p1].InterfaceList[j].interfacePartnerName.size();k++){
	       cout <<proList[p1].InterfaceList[j].interfacePartnerName[k]<<","<<proList[p1].InterfaceList[j].interfacePartnerIndex[k]<<'\t';
	  }//interface partners
	  cout <<endl;
      }//interfaces
      */
  }//proteins

  delete[] stoichvalue1;
  delete[] stoichvalue2;
  
}//end subroutine.

