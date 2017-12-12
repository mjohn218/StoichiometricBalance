#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double calc_frac_spec(int Nproteins, double *X, int *rxn_specific, int *prod_specific, int maxspecific, int **Rlist, double *X0, double &fracavg);
double calc_frac_spec_write(int Nproteins, double *X, int *rxn_specific, int *prod_specific, int maxspecific, int **Rlist, double *X0, double &fracavg, ofstream &outfile);
