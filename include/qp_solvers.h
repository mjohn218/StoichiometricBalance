#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

#include "QpGenData.h"
#include "QpGenVars.h"
#include "QpGenResiduals.h"
#include "GondzioSolver.h"
#include "QpGenSparseMa27.h"

using namespace std;

void qp_init(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *randvals);
void qp_solve(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *H, double *ZA, double *Q2, double ascale, int *p_home);
void qp_solve1(int Npro, int Nif, int Ncomplex, int &numstart, Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, constrainParms plist, int Nconstrain, int *constrain, double *abund, double * Q, double *c, double *xlow, double *xupp, char *ixlow, char *ixupp, int *irowQ, int *jcolQ, double *dQ, double *clow, double *cupp, char *iclow, char *icupp, int *irowC, int *jcolC, double *dC, double *randvals, double *H, double *ZA, double *Q2, double ascale, int *p_home, double avg_degree);
