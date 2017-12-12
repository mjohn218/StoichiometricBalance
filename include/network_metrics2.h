
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

double network2_metric_incom(int nwhole,  ppidata *ppi, int **Eref, int **Epred);
double network_metric(int nwhole, Protein *wholep, ppidata *ppi, int *numpartners, int **Speclist, Protein *refnet, int *refpartners, int **reflist, double alpha);
void locate_interfaces(int p1, int p2, Protein *wholep, int *numpartners, int **Speclist, int &p1if, int &p2if);
int find_islands(int N, int *Adj0, int **modlist, int *modsize);
