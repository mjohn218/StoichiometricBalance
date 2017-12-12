#include <fstream>
#include <iostream>
#include <string>

using namespace std;

double avg_rmsd_from_random_dists(int Npro, double *yeast_abund, double *cumul, double *abund, double *randvals, int Nconstrain, int *constrain, string *genid, int nruns, ofstream &rtotfile);
