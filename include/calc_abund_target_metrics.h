#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>

using namespace std;


double calc_msd_abund_target(int Npro, double *proconc, double *abund, ofstream &outfile, string *genid);
double calc_msd_nowrite_abund_target(int Npro, double *proconc, double *abund);
double calc_dist_entropy_totarget(int Npro, double *proconc, double *abund, string *genid);
double calc_dist_entropy_totarget_exclude_constrains(int Npro, double *proconc, double *abund, string *genid, int Nc, int *constrain );
double calc_corrcoef_totarget(int Npro, double *proconc, double *abund, string *genid, double &lgR);
double calc_corrcoef_totarget_exclude_constrains(int Npro, double *proconc, double *abund, string *genid, int Nc, int *constrain, double &lgR);
double calc_normeuclid(int Npro, double *proconc, double *abund);
double calc_mape(int Npro, double *proconc, double *abund);
double calc_chisquare(int Npro, vector<double>proconc, vector<double>abund);
double calc_dist_jensenshannon(int Npro, vector<double>proconc, vector<double>abund);
