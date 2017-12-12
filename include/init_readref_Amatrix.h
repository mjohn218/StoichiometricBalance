#include <fstream>
#include <iostream>

using namespace std;

void build_Amatrix_refProIIN(int Nif, int Nedge, double *A, int *e1num, int *e2num, Protein *wholep, int *numpartners, int **Speclist);
void init_AmatrixPPI(int Npro, int Nedge, int Nif, double *A, int *e1num, int *e2num, Protein *wholep, int *e1int, int *e2int);
