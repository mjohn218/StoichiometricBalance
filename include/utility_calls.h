#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <sys/time.h>
#include "rand_gsl.h"
#include <stdlib.h>

using namespace std;


#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define MI 7
#define NSTACK 50


double GaussV();
void indexx(unsigned long n, double arr[], unsigned long indx[]);
int trand();


