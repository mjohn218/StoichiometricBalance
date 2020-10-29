#include <vector>
#include <string>
#include "ProteinClass.h"
#include "pro_classes.h"
#include "construct_net.h"
using namespace std;

void reset_arrays(vector<string> &knockDownNames, vector<int> &proKDID, vector<ProteinClass> &proListKD, vector<double> &realKD, vector<double>&iFaceBalancedKD, vector<double> &proBalancedKD, vector<double> &stdBalancedKD)
{
  knockDownNames.clear();
  proKDID.clear();
  proListKD.clear();
  realKD.clear();
  iFaceBalancedKD.clear();
  proBalancedKD.clear();
  stdBalancedKD.clear();
  
}
