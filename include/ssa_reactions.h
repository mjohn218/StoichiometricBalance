
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>



using namespace std;


void mem_to_sol(int c1, gillFullmol *bases, gillComplex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *sol_prod);
void sol_to_mem(int c1, gillFullmol *bases, gillComplex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *mem_prod);
void mem_to_sol_skippip(int c1, gillFullmol *bases, gillComplex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *sol_prod);
void sol_to_mem_skippip(int c1, gillFullmol *bases, gillComplex *ind_com, double *X, int *ihome, int **myrxn, int **Rlist, int *p_home, int *mem_prod);
void break_complex_nocrds(int p1, int mu, int kind, gillFullmol *bases, int **Rlist, int *ihome, gillComplex *ind_com, Parms &plist, int p2, int i1, int i2, int *p_home, int **myrxn);
void break_complex_nocrds_skippip(int p1, int mu, int kind, gillFullmol *bases, int **Rlist, int *ihome, gillComplex *ind_com, Parms &plist, int p2, int i1, int i2, int *p_home, int **myrxn);
void break_complex_nocrds_3clath(int p1, int mu, int kind, gillFullmol *bases, int **Rlist, int *ihome, gillComplex *ind_com, Parms &plist, int p2, int i1, int i2, int *p_home, int **myrxn);
void associate_nocrds(int p1,int p2, int mu, int i1, int i2, gillFullmol *bases, int **Rlist, int *ihome, gillComplex *ind_com, Parms &plist);
int associate_nocrds_skippip(int p1,int p2, int mu, int i1, int i2, gillFullmol *bases, int **Rlist, int *ihome, gillComplex *ind_com, Parms &plist, int *p_home);
void set_status_conc(ifstream &startfile, gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, double *indivconc);
void set_status_conc_store(ifstream &startfile,  gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, double *indivconc, int *nfreestart, int **fliststart, int ** sliststart);
void reset_status_conc_store(gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, double *indivconc, int *nfreestart, int **fliststart, int ** sliststart);
void set_status_conc_skippip(ifstream &startfile, gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, double *indivconc);
void sum_state(int Nprotypes, gillProtein *wholep, int *Nmyrxn, int **myrxn, int **Rlist, int **Del, double *X, int *rxntype);
//int determine_which_complex_merge(int p1, int p2, gillComplex *ind_com, int c1, int c2, gillFullmol *bases, int **myrxn, int **Rlist, int *ihome, int *p_home, Parms &plist);
int determine_which_complex_double(int p1, int p2, gillComplex *ind_com, int c1, int c2, gillFullmol *bases, int **myrxn, int **Rlist, int *ihome, int *p_home, Parms &plist);
int determine_which_complex_skippip(int p1, int p2, gillComplex *ind_com, int c1, int c2, gillFullmol *bases, int **myrxn, int **Rlist, int *ihome, int *p_home, Parms &plist);

void update_diffusion(int c1, gillComplex *ind_com, gillFullmol *bases);
void update_Nlipids(int c1, int c2, gillComplex *ind_com, gillFullmol *bases);
