#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

double fitness_constrain(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_sqr(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_sqr2(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_concn(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_norm(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_norm_iface(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_nosq(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_nosq_iface(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund);
double fitness_constrain_parts(int nwhole, Protein *wholep, double *indivconc, int Nconstrain, int *constrain, double *net_abund, double &part1, double &part2);
