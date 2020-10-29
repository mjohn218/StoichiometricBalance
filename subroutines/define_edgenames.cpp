#include <string>
#include <vector>
#include "ProteinClass.h"
#include "pro_classes.h"
#include "construct_net.h"

/*Count up all the IIN edges, and give names to each edge based on interfaces participating
  This subrouting will cause the program to exit if it finds that the edges don't match the total Edges used
  to construct the arrays.

*/
void define_edgenames(int &Ncomplex, int Nif, int *numpartners, int **Speclist, string *edgenames, int nTotalEdges, string *iNames)
{
  int i, j;
  for(i=0;i<Nif;i++){
    for(j=0;j<numpartners[i];j++){
      if(Speclist[i][j]>=i){//Avoid counting an edge twice
	
	  edgenames[Ncomplex]=iNames[i];
	  edgenames[Ncomplex].append("::");
	  edgenames[Ncomplex].append(iNames[Speclist[i][j]]);
	  Ncomplex++;    
      }
    }
  }
  cout <<"-------------------"<<endl;
  cout << "IIN edges: " << Ncomplex << endl;
  if(Ncomplex != nTotalEdges){
      cerr <<" ERROR, EXITING: TOTAL DEFINED EDGES IN SPECLIST DOES NOT MATCH NUMBER DEFINED ABOVE "<<endl;
      cerr<<" Nc: "<<Ncomplex<<" doesn't match: "<<nTotalEdges<<endl;
      exit(1);
  }
  
}
