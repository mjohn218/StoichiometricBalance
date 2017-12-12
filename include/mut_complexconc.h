#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>


using namespace std;

int trand();
void mutate_comp(int ncomplex, int maxswap, int minnum, double *complexmut);
void mutate_glob(int ncomplex, int globswap, int minnum, double *complexmut);
void mutate_single(int ncomplex, int maxswap, int minnum, double *complexmut);
void mutate_single_gauss(int ncomplex, int maxswap, int minnum, double *complexmut);
