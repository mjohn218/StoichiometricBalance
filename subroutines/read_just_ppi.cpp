#include "pro_classes.h"
#include "constrainParms.h"
#include "ProteinClass.h"
#include "read_proinput.h"

int read_just_ppi(int nwhole, ifstream &ppifile, ppidata *ppi)
{
  int i, j;
  int ig;
  int npart;
  int nedge=0;
  for(i=0;i<nwhole;i++){
    ppifile>>ig; 
    ppifile>> npart;
    ppi[i].nppartner=npart;
    nedge+=npart;
    for(j=0;j<npart;j++)
      ppifile >>ppi[i].pplist[j];
  }
  nedge/=2;
  return nedge;
}
