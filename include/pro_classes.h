
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>

#define MAXIFACE 200

#define MAXP 200


using namespace std;


#define EDIM 200
#define NBINS 16000

class Protein
{
public:
  int ninterface;
  int valiface[MAXIFACE];
};
class ppidata
{
public:
  int nppartner;
  int pplist[MAXP];
};

