#include <fstream>
#include <iostream>

using namespace std;

void write_msd_abund_target(int Npro, double *proconc, double *abund, ofstream &outfile, string *genid, double *randvals);
void write_copt_three(int Nedge, double *complexopt, double *restart, double *restart2, ofstream &outfile, int Nc, int *constrain );
void write_copt_pair(int Nedge, double *complexopt, double *restart, ofstream &outfile, int Nc, int *constrain );
void write_copt(int Nedge, double *complexopt, ofstream &outfile, int Nc, int *constrain);
void write_prot(int nwhole, Protein *wholep, double *conc, ofstream &outfile, double &nprot);
void write_complex(int Nedge, double *conc, int *e1num, int *e2num, int *e1int, int *e2int, ofstream &outfile);
void write_conc(int nmol, double *conc, ofstream &outfile);
void write_parms(constrainParms &plist);
