#include "ProteinClass.h"
#include "constrainParms.h"

void read_proteins_from_file(int &Npro, ifstream &expfile, constrainParms &plist, vector<string> &genid, vector<double> &abund,  vector<ProteinClass> &proList);
