#include <vector>
#include "pro_classes.h"
#include "constrainParms.h"
#include "ProteinClass.h"
#include "read_proinput.h"

void remove_duplicate_edges(int &NedgeRead, vector<string> &edge1Interface, vector<string> &edge2Interface, vector<string> &edge1Protein, vector<string> &edge2Protein, vector<string> &stoichiometry)
{
 
  int i;
  bool duplicate=false;
  for(i=1;i<NedgeRead;i++){
      //look at all previous edges for duplicates
      for(int prev=0;prev<i;prev++){
	  duplicate=false;
	  if(edge1Protein[prev]==edge1Protein[i] && edge2Protein[prev]==edge2Protein[i] && edge1Interface[prev]==edge1Interface[i] && edge2Interface[prev]==edge2Interface[i])
	      {
		  cout <<"READ IN THE SAME EXACT INTERACTION TWICE, DELETING ONE OF THEM. Lines: "<<prev+2<<" and "<<i+2<<endl;
		  /*Replace the previous edge with the last edge in the list, then decrease the size of the array.*/
		  int last=NedgeRead-1;
		  edge1Protein[prev]=edge1Protein[last];
		  edge2Protein[prev]=edge2Protein[last];
		  edge1Interface[prev]=edge1Interface[last];
		  edge2Interface[prev]=edge2Interface[last];
		  stoichiometry[prev]=stoichiometry[last];
		  NedgeRead--;
		  duplicate=true;

	      }
	  if(duplicate==false){
	      //check the reversed A B interaction, but not if B A is same for this pair.
	      if(edge1Protein[prev]==edge2Protein[i] && edge2Protein[prev]==edge1Protein[i] && edge1Interface[prev]==edge2Interface[i] && edge2Interface[prev]==edge1Interface[i])
		  {
		      cout <<"READ IN THE SAME EXACT INTERACTION TWICE, (REVERSED A AND B COLUMNS). DELETING ONE OF THEM. Lines: "<<prev+2<<" and " <<i+2<<endl;
		      /*Replace the previous edge with the last edge in the list, then decrease the size of the array.*/
		      int last=NedgeRead-1;
		      edge1Protein[prev]=edge1Protein[last];
		      edge2Protein[prev]=edge2Protein[last];
		      edge1Interface[prev]=edge1Interface[last];
		      edge2Interface[prev]=edge2Interface[last];
		      stoichiometry[prev]=stoichiometry[last];
		      NedgeRead--;
		      
		  }
	  }
      }
  }
  cout <<" -------------------- "<<endl;
  cout <<"AFTER CHECKING FOR DUPLICATES, EDGES READ IN IS NOW: "<<NedgeRead<<endl;
  cout <<" -------------------- "<<endl;

}
