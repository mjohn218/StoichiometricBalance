#include "pro_classes.h"
#include "constrainParms.h"
#include "ProteinClass.h"
#include "read_proinput.h"

void read_short_parms(ifstream &parmfile, constrainParms &plist)
{
  
  parmfile >>plist.min_complex;
  parmfile.ignore(400,'\n');
  parmfile >>plist.ascale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.nruns_rand;//If this is zero, it will still evaluate balance for REAL copy numbers
  parmfile.ignore(400,'\n');
  parmfile >>plist.flagread;//this is only relevant for randomized copy number runs.
  parmfile.ignore(400,'\n');
  parmfile >>plist.max_complex;//If this is zero, will allow max up to infinity. Should at least be big enough to allow for randomized copy numbers to span the full maximal distribution. 
  parmfile.ignore(400,'\n');
  
}
