#ifndef __CLASS_DEFS_H
#define __CLASS_DEFS_H

#include <fstream>
#include <string>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>

#define MAXIFACE 20
#define MAXPRTNER 30
#define MAXCOMPLEX 50000
#define MAXRXN 200
#define MAXOVERLAP 40


using namespace std;

/* class Mbind */
/* { */
/* public: */
/*   int p1start; */
/*   int p2start; */
/*   int ncol; */
/*   int nrow; */
/* }; */
class gillProtein
{
public:
  int ninterface;
  int valiface[MAXIFACE];
  int npropart;
  int propart[MAXPRTNER];
  double Dx;
  double Dy;
  double Dz;
  
  
};
class gillFullmol
{
public:
  
  int protype;
  int ninterface;
  int istatus[MAXIFACE];
  int mycomplex;
  int Nlipid;
  int nbnd;
  int nfree;
  int freelist[MAXIFACE];
  int bndlist[MAXIFACE];
  
  int partner[MAXIFACE];
  double Dx;
  double Dy;
  double Dz;
    
};
class gillComplex
{
public:
  int mysize;
  int Nlipid;
  double Dx;
  double Dy;
  double Dz;
  int plist[MAXCOMPLEX];
  

};
class Parms
{
public:
  int Nprotypes;
  int Nifaces;
  double Nit;
  int restart;
  int statwrite;
  int configwrite;
  int grwrite;
  double deltat;
  double V;
  double X0total;
  int Nspecies;
  int Nrxn;
  double mass;
  double D;
  int nspec_complex;
  double maxsep2;
  int ntotalcomplex;
  double xboxl;
  double yboxl;
  double zboxl;
  int Ntotalmol;
  int Natom;
  int Natomwrite;
  double pretrans;
  double prerot;
  int pclath;
  int nloop;
};
class ProductName
{
 public:
  int index;
  int i1;
  int i2;
  string cname;
  string sname;
  string mname;
  string ifacename;
  string prodname;
};
#endif
