#include <iostream>
#include <fstream>

using namespace std;

double opt_complexconc_mc(constrainParms plist, int Nif, int Npro, int Ncomplex, Protein *wholep, double *A, double *indivconc, double *complexconc, double *complexopt,int *constrain, double *abund,  double *yeast_abund, double *cumul, int Nconstrain, int mtypes, double *moveprob, string *genid, ofstream &idfile, ofstream &matchfile, ofstream &compfile, int nc, int nr, Metrics &mets, Metrics &rndmets);
double opt_complexconc_mc_nowrite(constrainParms plist, int Nif, int Npro, int Ncomplex, Protein *wholep, double *A, double *indivconc, double *complexconc, double *complexopt,int *constrain, double *abund,  double *yeast_abund, double *cumul, int Nconstrain, int mtypes, double *moveprob, string *genid, Metrics &mets, Metrics &rndmets);
