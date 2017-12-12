#include "pro_classes.h"
#include "constrainParms.h"
#include "write_outs.h"

void write_parms(constrainParms &plist)
{
  cout <<"parameters read in.. "<<endl; 
  cout <<"N whole proteins: "<<plist.nwhole<<endl;
  cout <<"Min Number of molecules to allow: "<<plist.min_complex<<endl;
  cout <<"Nedges: "<<plist.Nedge<<endl;
  cout <<"Nconstrain: "<<plist.Nconstrain<<endl;
  cout <<"ascale: "<<plist.ascale<<endl;  
  
}
