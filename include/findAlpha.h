#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdlib> // generates randoom number
#include <ctime>   //time seed!
#include <numeric>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

int index (int row , int column, int width);
void printarray( int *array, int length);
void printarraydoub( double *array, int length);
void print2darray (int *TwoDarray, int rows, int columns);
vector<double> divideArray(vector<double> array, double divisor, int length);
double sum_array(vector<double> a, int num_elements);
int sum_arrayint(int *a, int num_elements);
vector<double> Cum_Sum(vector<double> array, int length);
double *ScaleFree(int Nodes, double alpha);
int CountOrphans(int *array, int length);
int * ZeroArray( int * array, int length);
vector<double> ZeroArrayDouble( vector<double>  array, int length);
vector<int>  locateIndexOfZeros(vector<int> Vect, int * array, int length);
void printVect(vector<int> vecArray);
int *NeighborArray( int *TwoDArray, int Nodes);
int Pick(double *ProbVect, int Nodes);
int leaveout(double *NormProbs, int Nodes, int Edges, double Alpha);
int *GenerateAlphaNetwork(int Nodes, int Edges, double alpha);
double findvals(int * array, int value, int length);
vector<double> PDFcalculate (int * Network, int dimentions);
vector<double> addTwoArrays(vector<double> array1, vector<double> array2, int length);
vector<double> AvgPDF (int Nodes, int Edges, double alpha, int numIters);
double chierror(vector<double> pmf1 , vector<double> pmf2, int length);
int indexofSmallestElement(vector<double> array, int size);
double findAlpha(double *pmfCompare , int Edges, int Nodes, int nwhole, int maxni, int PPIedge, int maxne, double **sampPMFs);
double findAlpha2(double *pmfCompare, int Edges, int Nodes);
void readPMFfile(ifstream &pmf_file, double **PMFs, int lines);
