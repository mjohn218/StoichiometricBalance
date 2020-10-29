#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>


using namespace std;


void read_ref_ProIIN(int nwhole, Protein *wholep, int &nrefinterface, ifstream &reffile, int *npartners, int **Speclist, int *p_home);
void read_cnstrn_parms(ifstream &parmfile, constrainParms &plist);
void read_balance_parms(ifstream &parmfile, constrainParms &plist);
void read_short_parms(ifstream &parmfile, constrainParms &plist);
void read_protlist(int nwhole, Protein *wholep, int nmol, int *p_home, ifstream &protfile);
int read_just_ppi(int nwhole, ifstream &ppifile, ppidata *ppi);
void read_edges_from_file(ifstream &iinfile, vector<string> &edge1Interface, vector<string> &edge2Interface, vector<string> &edge1Protein, vector<string> &edge2Protein, vector<string> &stoichiometry, vector<string> filenames);
void remove_duplicate_edges(int &NedgeRead, vector<string> &edge1Interface, vector<string> &edge2Interface, vector<string> &edge1Protein, vector<string> &edge2Protein, vector<string> &stoichiometry);
void read_copynumber_dist(char *copyNumberFile, vector<double> &bin, vector<double> &cdf, constrainParms &plist);
void read_knockdowns(ifstream &kdfile, vector<string> &knockDowns, int Npro, vector<ProteinClass> proList, vector<int> &proKDID, int &nKD);
