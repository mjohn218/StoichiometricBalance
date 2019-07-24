#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;


void read_ref_ProIIN(int nwhole, Protein *wholep, int &nrefinterface, ifstream &reffile, int *npartners, int **Speclist, int *p_home);
void read_cnstrn_parms(ifstream &parmfile, constrainParms &plist);
void read_balance_parms(ifstream &parmfile, constrainParms &plist);
void read_protlist(int nwhole, Protein *wholep, int nmol, int *p_home, ifstream &protfile);
int read_just_ppi(int nwhole, ifstream &ppifile, ppidata *ppi);
