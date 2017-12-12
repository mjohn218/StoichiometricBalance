#include <iostream>
#include <fstream>
using namespace std;

void write_conc_vs_time(ofstream &corrfile1, ofstream &corrfile2, int it, double curr_time, double *X, int Nifaces, int Nspecies);
void write_conc_vs_time_nonzero(ofstream &corrfile1, ofstream &corrfile2, int it, double curr_time, double *X, int Nifaces, int Nspecies, ProductName *complexnames);
