#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>

#define MAXPRTNER 2000 

using namespace std;

class nParms{
public:
  int Nstep;
  double beta;
  double T1;
  int pcut;
  int twrite;
  int flagself;
};
void mutate_net_degree_preserve(int N, int *npartner, int **Speclist, int *Adj, int &chg );
void accept_net(int N, double newfit, double &oldfit, int *Adj, int *Tempadj, int **Speclist, int *temppartner, double &globopt, int *Globadj);
void gen_and_opt_net(int N,int Totedge, double alpha, int *Adj, int **Speclist, int *numpartner, double *histavg, double &cglobal, double &clocal, double &distmean, double &diststd, double &cglobalopt, double &clocalopt, double &diststdopt, nParms &plist, double &normavg, double &normstd, double &normavgo, double &normstdo, double *histfour, double *histfouropt,  int nr, double &ctri, double &ctrio, int flagself);

double nbor_fitness(int N, int *Adj, int **Speclist , int pcut, double beta);
double ratio_fitness(int N, int *Adj, int **Speclist , int pcut, double beta);
void optimize_network_conn(int N, nParms &plist, int *Adj, int **Speclist, int *numpartner, int nr);
void starting_network(int N, int Totedge, double alpha, int *Adj, int **Speclist, int *numpartner, double *hist, int nr, int *self, int &Nself);
void starting_network_noself(int N, int Totedge, double alpha, int *Adj, int **Speclist, int *numpartner, double *hist, int nr);


void print_network_edge_list(int N, int *Adj, char *filename);
void print_network_node_then_partners(int N, int *numpartners, int **Speclist, int nr);
void print_network_node_then_partners(int N, int *Adj, int nr);
