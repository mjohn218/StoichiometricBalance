#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>


using namespace std;


class rdParms
{
public:
  int nmol;
  int Nsite;
  double V;
  double conc_start;
  int readconc;
  double Maxtime;
  int checkwrite;
  int stepwrite;
  double Kdnon;
  double Kdspec;
};


void read_parms_net(ifstream &parmfile, rdParms &plist);

double gillespie_allrxn(int Nrxn, int nmol, int maxspecific, double *X, int **Rlist, int*rxntype, double *c, double *a, int **Del, int **myrxn, int *prod_specific, int *rxn_specific, int Nmyrxn, double MaxTime, int *Npart, int time, int step, GillParms &plist, double *cstart);

void set_rxntype(int Nrxn, int *Npart, int **Del, int *rtype);
void set_binary_reactions(int Nproteins, int *Npart, int **Rlist, int **Del);
void write_parms(rdParms &plist);
void set_initial_conc(int Nproteins, int Nspecies, double *X0, double *X, int numstart);
void write_final(int Nspecies, double *X, ofstream &ffile);
void read_initial(int Nspecies, double *X, ifstream &concfile);
void set_all_reactions(int Nproteins, int *Npart, int **Rlist, int **Del, int **myrxn);


void match_reactions(int Nproteins, int **Rlist, int **Del, int *Npart, int *numpartner, int **Speclist, int *rxn_specific, int *prod_specific, int &maxspecific);
void write_spec(int Nproteins, double *X, int *rxn_specific, int *prod_specific, int maxspecific, ofstream &ffile, int **Rlist);
void write_misbind(int Nproteins, double *X, int **Rlist, int Nspecies, ofstream &ffile);

void initialize_rates_nonspecific(int Nproteins, double *k, double Kdnon);
void read_initial_conc(int Nspecies, double *X, ifstream &concfile, double V);


/*assigning rates for gillespie to use from matrices of rates*/
void get_all_rates(int Nproteins, double *k, double *k_seq);
void get_rates_misbind(int Nproteins, double *k, double **kmis_seq, int Totmisbind);
void set_misbind_reactions(int Nproteins, int *Npart, int **Rlist, int **Del, int Totmisbind, int **myrxn);
void write_rate(ofstream &rfile, int nmol, double *k_seq, int flagmis, double **kmis_seq, int Totalmisbind);
