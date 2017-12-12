
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <stdlib.h>



using namespace std;

class GillParms
{
public:
  int Vsphere;
  int nwhole;
  int Nif;
  double Maxtime;
  double V;
  double X0total;
  int Nrxn;
  int Nspecies;
  int stepwrite;
  int checkwrite;
  int Nit;
  double delt;
  int Nruns;
  int Ntotalmol;
  int ntotalcomplex;
  int Nprotypes;
  int Nifaces;
  /*more parms for dealing with sequence binding energies*/
  int flagmis;
  int flagrot;
  int Nsite;
  double Kdspec;
  double kbt;
  double mincom;
};

