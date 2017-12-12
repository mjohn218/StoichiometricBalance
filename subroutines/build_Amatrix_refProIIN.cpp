#include "pro_classes.h"
#include "init_readref_Amatrix.h"

void build_Amatrix_refProIIN(int Nif, int Nedge, double *A, int *e1num, int *e2num, Protein *wholep, int *numpartners, int **Speclist)
{
  int i, j;
  int part;
  int t=0;
  
  int nc=0;
  int p1, p2;
  int i1, i2;
  int k, s;
  for(i=0;i<Nedge;i++){
    p1=e1num[i];
    p2=e2num[i];
    //find the interfaces for this edge
    
    //    cout <<"p1: "<<p1<<" p2: "<<p2<<endl;
    for(s=0;s<wholep[p1].ninterface;s++){
      //loop over the interfaces on the first protein
      i1=wholep[p1].valiface[s];
      for(j=0;j<wholep[p2].ninterface;j++){
	//loop over the interfaces on the second protein
	i2=wholep[p2].valiface[j];
	/*See if these 2 interfaces interact*/
	for(k=0;k<numpartners[i1];k++){
	  //the list of i1's partners
	  if(Speclist[i1][k]==i2){
	    //then this is the edge, break out of the loop
	    //	    cout <<"match: "<<i1<<' '<<i2<<endl;
	    k=numpartners[i1];
	    j=wholep[p2].ninterface;
	    s=wholep[p1].ninterface;
	    
	  }	      
	}
      }
    }
    
    
    
    
    //cout <<"protein 1: "<<p1<<" myiface: "<<i1<<endl;
    //cout <<"protein 2: "<<p2<<" myiface: "<<i2<<endl;
	
    if(i1==i2){
    
      A[i*Nif+i1]=2;
      
    }
    else{
      
      //in this case each interface has a single partner
      //rows are for each interface, columns are for each complex
      A[i*Nif+i1]=1;
      A[i*Nif+i2]=1;
    }
    nc++;
    //cout <<"ncomplex: "<<nc<<" p1: "<<p1<<" p2: "<<p2<<endl;
    
  }
   

}
