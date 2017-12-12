#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdexcept>

using namespace std;

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define MI 7
#define NSTACK 50

void four_motif(int N, int *Adj, double **fourmer, double *histfour, double **hist, int *typelist);
void cluster_cof(int N, int &open, int &closed, double &avg, int  **Speclist, int *numpartner);
void cluster_cof_self(int N, int &open, int &closed, double &avg, int  **Speclist, int *numpartner, int *self, double &selfavg);
void neighbor_conn(int N, int *numpartner, int **Speclist, double  *conn, double *ctot, double &totavg, double &std, double &normavg, double &normstd);
void indexx(unsigned long n, double arr[], unsigned long indx[]);
void order_nodes(int Nnode, int *origlist, int *ordered);
void three_ways(int i, int j, int k, int l, int *Adj, int N, int &t, double **fourmer, double **hist, int *origlist, int *typelist);
int check_repeat(int Nmer, double **fourmer);
void find_type(int i, int j, int k, int l, int *Adj, int N, double **hist, int tak);
void grid_cof(double *gcs, int N, int *numpartners, int **Speclist);
double cluster_cof_mod(int i1, int N, int *numpartners, int **Speclist);
