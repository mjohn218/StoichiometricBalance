
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string>
#include "write_conc_vs_time.h"

using namespace std;

void check_reactions(Parms plist, int *numpartners, int **Speclist, int *Nmyrxn, int **Rlist, int **myrxn, int Nspecies);

void read_reactions_complex(ifstream &rxnfile, int Nrxn, int **Rlist, int *Nmyrxn, int **myrxn, int *rxntype, double *k, int Nspecies, int **Del, int *Npart, int *Ncoup, int **mycoupled, int *rxn2D, int Nifaces, int *mem_prod, int *sol_prod, ProductName *complexnames);

/*below is with also RD states*/
void read_protlist_twostate(int nwhole, gillProtein *wholep, int Niface, int *p_home, ifstream &protfile, int *ihome, string *pronames);
void read_network_comment(Parms &plist, int *numpartner, int **Speclist, ifstream &infile);
