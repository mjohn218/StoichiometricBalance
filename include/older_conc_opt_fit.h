#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

void conc_fitness(int nwhole, int ninterface, Parms &plist, int *numpartners, int **Speclist, Protein *wholep, int *p_home, double *concfit, double *Amat, double *indivconc, double *complexconc, double *complexmut, double *complexopt);
void optimize_complex(Parms &plist, int ncomplex, Protein *wholep, int *p_home, double *indivconc, double *A, double *complexconc, double *complexmut, double *complexopt, int &npfinal, double &fvec1);
double fitness(int nwhole, Protein *wholep, double *indivconc);
