#ifndef __GILL_SUBS_H
#define __GILL_SUBS_H

#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>



using namespace std;



void read_Gillparms(ifstream &parmfile, GillParms &plist);
void read_Gillparms_complex(ifstream &parmfile, GillParms &plist);

double hval(int rxn, int **Rlist, int *rtype, double *X);
void set_cvals(int Nrxn, double *k, double *c, int *rtype, double V);
void set_cvals_mem(int Nrxn, double *k, double *c, int *rtype, double Vmol, double V2D, int *rxn2D);
void update_cvals_rho(int Nrxn, double *k, double *c, int *rtype, double Vmol, double V2D, int *rxn2D, double D, double sigma, double ka, double *X, int **Rlist, double kb);
void update_cvals_rho_vec(int Nrxn, double *k, double *c, int *rtype, double Vmol, double V2D, int *rxn2D, double *Dvec, double *sigma,  double *X, int **Rlist, double *Kdiss);

double gillespie(int Nrxn, int nmol, double *X, int **Rlist, int *rxntype, double *c, double *a, int **Del, int **myrxn, int *Nmyrxn, double MaxTime,  int *Npart, int time, int step, GillParms &plist, double *cstart );
double gillespie_avg(int Nrxn, int nmol, double *X, int **Rlist, int*rxntype, double *c, double *a, int **Del, int **myrxn, int *Nmyrxn, double MaxTime,  int *Npart, int time, int step, GillParms &plist, double *cstart, double *Aavg, double *Ntimes, double delt);

void read_reactions(ifstream &rxnfile, int Nrxn, int **Rlist, int *Nmyrxn, int **myrxn, int *rxntype, double *k, int Nspecies, int **Del, int *Npart);
void read_startlist(ifstream &scfile, double *concstart, double *indivstart, double V, int nmol);
void read_molnumber(ifstream &scfile, double *concstart, double *indivstart, double V, int nmol);
void read_net(int Nprot, int *npartners, int **neighbor, ifstream &infile);




#endif
