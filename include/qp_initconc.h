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

void qp_initconc(int Npro, int Nif, int Ncomplex,  Protein *wholep, double *yeast_abund, double *cumul, double *indivconc, double *complexconc, double *A, double avg_degree);
