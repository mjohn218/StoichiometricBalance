#include <vector>
#include <string>
#include "constrainParms.h"
#include "ProteinClass.h"
#include "pro_classes.h"
#include "read_proinput.h"


/*If we want to calculate a p-value (OPTIONAL) of the observed stoichiometry relative to randomized copy numbers, 
  need to read in a distribution of copy-numbers from the corresponding cell-type.
  
*/
void read_copynumber_dist(char *copyNumberFile, vector<double> &bin, vector<double> &cdf, constrainParms &plist)
{
  int t=1;
  ifstream distfile(copyNumberFile);
  double btmp;
  double ctmp;
  bin.push_back(0);//first value of copy number bins
  cdf.push_back(0);//first value of CDF

  /*Read in conc distribution*/
  cout <<"---------------"<<endl;
  cout <<"READING IN CELL TYPE COPY NUMBER DISTRIBUTION FROM THE FILE: "<<copyNumberFile<<endl;  
  if(plist.flagread==1)
    cout<<"For P-VALUE: Using Copy number sampled from the distribution in: "<<copyNumberFile<<endl;
  else
    cout<<"For P-VALUE: Using Copy Numbers shuffled from the real COPY NUMBERS in Protein file: "<<endl;
  distfile.ignore(400,'\n');//ignore column name titles.
  while(distfile >> btmp >> ctmp){
    bin.push_back(btmp);
    cdf.push_back(ctmp);
    t++;
  }
  distfile.close();
  cout <<"TOTAL LINES in Copy Number Distribution for Cell Type: "<<t<<endl;

  /*make sure max_complex is big enough to accomodate the large protein copy numbers in this distribution.*/
  if(plist.max_complex < bin[bin.size()-1] && plist.max_complex>0){
    cout<<"---------------------WARNING-----------------------"<<endl;
    cout <<" WARNING: max_complex value is too low--limits the complexes and they might not be able to sum up to maximum copy numbers observed in this cell type. "<<'\t';
    cout <<"setting max_complex to the largest copy number in the cell-type: ";
    plist.max_complex=bin[bin.size()-1];
    cout<<plist.max_complex<<endl;
    cout<<"--------------------------------------------"<<endl;
  }
  
}
