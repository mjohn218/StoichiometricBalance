#include "pro_classes.h"
#include "constrainParms.h"
#include "read_proinput.h"

void read_cnstrn_parms(ifstream &parmfile, constrainParms &plist)
{
  parmfile >>plist.nwhole;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nedge;
  parmfile.ignore(400,'\n');
  parmfile >>plist.min_complex;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nconstrain;
  parmfile.ignore(400,'\n');
  parmfile >>plist.flagread;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Ncsets;
  parmfile.ignore(400,'\n');
  parmfile >>plist.nruns_rand;
  parmfile.ignore(400,'\n');
  parmfile >>plist.ascale;
  parmfile.ignore(400,'\n');
  parmfile >>plist.nstepMC;
  parmfile.ignore(400,'\n');
  parmfile >>plist.mut_change;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Tone;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nt;
  parmfile.ignore(400,'\n');
  parmfile >>plist.scaleT;
  parmfile.ignore(400,'\n');
  parmfile >>plist.Nruns;
  parmfile.ignore(400,'\n');
  

  
}
