#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;


void build_Amatrix(int nmol, double *A, int * numpartners, int **Speclist);
void build_Amatrix_with_Stoich(int nmol, double *A, int * numpartners, int **Speclist, double **Stoichlist);
