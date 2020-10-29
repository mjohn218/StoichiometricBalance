#include <string>
#include <sstream>
#include <vector>

using namespace std;

void sample_dist(int npro, vector<double> bin, vector<double> cdf, vector<double> &rval);
void shuffle_abund(int npro, vector<double> abund, vector<double> &rval, int Nconstrain, vector<int> constrain);
void get_variances(int Npro, int Nif, Protein *wholep, vector<double> R, double *ivars, double *imean);
