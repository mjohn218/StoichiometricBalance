/*

./stoich_balance.exe parms.inp prolist_abund_constrain.inp edgelist.inp knockDownList  >std_out.out
(file is optional)


July 2019, updated input/output formats.

Inputs:
- parms.inp                set of parameters
- prolist_abund_constrain.inp          Protein Names, abundance of proteins, isConstrained?
- edgelist.inp             list of protein1, interface1, protein2, interface2 interactions.
OPTIONAL: This is only needed if you want to calculate the p-value for deciding if the observed copy numbers are balanced, with significance.
- CellType_AllProteinCopyNumber_Distribution.txt    Distribution of protein copy numbers for this cell type. Used only to be used to sample random copy numbers for p-value calculation. Two columns: copy numbers, Cumulative distribution function.

Outputs:
1- "RealCopyNumbers_vsBalancedCopyNumbers.out"  The balanced copy number solution, given the real, observed copy numbers.
2- ChiSquaredDistance.txt   The CSD between the starting and balanced copy numbers, for the Real copy numbers and those initialized based on the CellTypeDistirbution.
3- JensenShannonDistance.txt The JSD between the starting and balanced copy numbers, for the Real copy numbers and those initialized based on the CellTypeDistirbution.
4- Vars.out  The variances of the interface copy numbers across each protein, for each solution. (Real + all randomized).
5- Means.out The mean of copy numbers for each protein, for each solution (Real + all randomized).
6- stdout.out  This tells you what parameters you inputted for the run, the IIN that was constructed based on the proteins listed in the prolist_abund_constrin.inp file, protein properties, number of constraints on copy numbers, and success of each balancing run.



OPERATION:

Code minimizes a quadratic function of the form

min_x (Ax-C0)'Z(Ax-C0) + (Ax)'H(Ax)

See manuscript Methods section for description of the matrices

The first C0 used is the real copy numbers, defined in node_abund.inp. Only the copy numbers of constrained proteins are used; for unconstrained proteins Z(i,i) = 0 and thus the term does not add to the equation. C_balanced is then calculated via

A*x_min = C_balanced

The chi-square distance and Jensen-Shannon distance between C0 and C_balanced is calculated.

Next, a random set of copy numbers for C0 is sampled Nit times, and the process is repeated each time. The CSDs and JSDs for each set of random copy numbers is printed out to be compared to the real C0. This code does not calculate the p-value, however. The copy numbers can be shuffled instead of randomly sampled by setting the parameter [NAME] in parms.inp to [VALUE].

Also calculated is the distance from C0 to the interfaces on P_balanced; that is, when the interface copy numbers on the same protein have been averaged. E.g. if in C_balanced the copy numbers of three interfaces on the same protein are [100 100 130], then for P_balanced they will be [110 110 110]. For ease, in the code these vectors are referred to as:

R:        C0
Rp:       C_balanced
Rpp:      P_balanced

Also outputed is the balanced solution for the real C0, so you can read the balanced protein copy numbers and compare them to the real copy numbers.

*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <sstream>
#include <vector>

#include "md_timer.h"
#include "pro_classes.h"
#include "read_protein_class.h"
#include "metrics_class.h"
#include "matmultiply.h"
#include "qp_solvers.h"
#include "functionCalls.h"
#include "initialize_qp.h"
#include "solve.h"

#include "read_proinput.h"
#include "construct_net.h"
#include "read_yeast_conc.h"
#include "write_outs.h"
#include "fill_iface_conc.h"
#include "gen_fullrand_dist.h"
#include "fitness_constrain_to_target.h"

#include "avg_rmsd_from_random_dists.h"
#include "build_Amatrix.h"

#include "calc_abund_target_metrics.h"
#include "gen_constraints.h"
#include "gen_rand_dist.h"
#include "even_protein.h"
#include "init_iin_net.h"

using namespace std;

struct MD_Timer totaltime;
struct MD_Timer inittime;

 
int main(int argc, char *argv[])
{
  
  /*Define the reactants, and all the species involved in the reactions
   *this includes all the products, including the misbinding products
   */

  int i, j;
  timeval tim;
  gettimeofday(&tim, 0);
  double t1=tim.tv_sec+tim.tv_usec;
      
  int seed=t1;
  srand(seed);
  
  ifstream parmfile(argv[1]);
  constrainParms plist;
  plist.max_complex=0;//initialize in case it is not read in.
  vector<string> filenames;
  filenames.push_back(argv[1]);//parmfile
  filenames.push_back(argv[2]);//protein list
  filenames.push_back(argv[3]);//edgelist
  
  /*READ IN PARAMETERS FROM ARGV[1] FILE*/
  read_short_parms(parmfile, plist);
  write_parms(plist);
    

  int Npro=1;

  vector<string> genid;
  vector<double> abund;

  
  vector<ProteinClass> proList;
  

  string text;
  int t=0;
  
  /*READ IN LIST OF PROTEINS AND KNOWN CONCENTRATIONS */
  ifstream expfile(argv[2]);
  read_proteins_from_file( Npro, expfile,  plist,  genid,  abund,  proList);
  plist.nwhole=Npro;
  
  /*READ IN NETWORK OF INTERACTIONS BETWEEN PROTEIN INTERFACES */
  ifstream iinfile(argv[3]);
  /*Read P1 I1 P2 I2*/
  vector<string> edge1Interface;
  vector<string> edge2Interface;
  vector<string> edge1Protein;
  vector<string> edge2Protein;
  vector<string> stoichiometry;
  
  read_edges_from_file(iinfile, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames);
  int NedgeRead=edge1Interface.size();//this can be reduced if duplicates are found.
  
  /*REMOVE DUPLICATES FROM THE INTERFACE INTERACTION NETWORK, THESE CAN OCCUR IF YOU READ PAST THE END OF FILE.*/
  remove_duplicate_edges(NedgeRead, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry);

  /*READ IN A LIST OF PROTEINS TO KNOCK DOWN*/
  ifstream knockDownFile;//(argv[4]);
  int nProKD=0;
  vector<int>proKDID;
  vector<string> knockDownNames;
  
  vector<ProteinClass> proListKD;
  int NproWithKD;//=Npro-nProKD;
  //proListKD.reserve(NproWithKD);

  bool loopAll=false;
  if(strncmp(argv[4], "-all",4)==0){
    cout<<" Perform knockdown of each individual protein "<<endl;
    loopAll=true;
  } else{
    knockDownFile.open(argv[4]);
    read_knockdowns(knockDownFile, knockDownNames, Npro, proList, proKDID, nProKD);
    create_new_prolist(proList, proListKD,  proKDID,  nProKD);
  }
  /*THIS CODE WILL NOT CALCULTE P-VALUES, SO IF NRUNS IS SET, PUSH TO ZERO.*/
  cout<<" NO CELL TYPE DISTRIBUTION FILE READ IN. NOT PERFORMING STOICHIOMETRIC BALANCE FOR ANY RANDOMIZED COPY NUMBERS--NO P-VALUE CALCULATION."<<endl;
  
  /*HERE, WE CAN PERFORM THE WHOLE ANALYSIS BASED ON THE INPUT FILES.
    FOR COMPARISONS, PRE AND POST-KNOCKDOWN, WE WILL NEED:
    R, R', R''. 
    The names of each R, which are stored in proList.
    With those distributions, we can normalize them, and calculate distances
    and changes per each protein. 
    So for each knockdown, we would like a distinct set of R, R' and R'' values, 
    and a distinct set of proList generated. 
  */
  ofstream distanceFile;
  distanceFile.open("ChiSquared_JensenShannonDistances.out", ios::out | ios::trunc);

  string label ="Original";
  distanceFile<<label << " Real to R'\tReal to R''\tR' to R''\n";
  string label2 = "KD";
  for(int i=0;i<nProKD;i++)
    label2=label2+":"+knockDownNames[i];
  string fappend=label2;
  if(loopAll==true)
    fappend="loopAll";
  ofstream deltaFile;
  deltaFile.precision(12);
  char fname[300];
  sprintf(fname, "DeltaOriginal_vsKDsolutions_%s.txt",fappend.c_str());//,Npro, label.c_str());//label is a string.
  deltaFile.open(fname);
  deltaFile << "KnockedDown\t"<<"CSDKDvsREALNif\t"<<"JSDKDvsREALNif\t"<<"CSDOrigVsKDifaceNif\t"<<"CSDOrigVsKDProNif\t"<<"JSDOrigVsKDIfaceNif\t"<<"JSDOrigVsKDProNif\t"<<"ChiRBratioNPRO\t"<<"JSDRBratioNPRO"<<endl;
  vector<double> Real;
  vector<double> iFaceBalanced;//Balanced solution--interfaces on a protein are not necessarily the same conc.
  vector<double> proBalanced; //Balanced solution after each protein is assigned an average over all interfaces
  vector<double> stdBalanced; //stand. deviation within a protein of interface copies.
  vector<double> distMetrics; //return the CSD and JSD for the balanced vs real copies.
  
  /*Below returns altered Real, IfaceBalanced, proBalanced, varBalanced, proList*/
  solve_balanced_from_inputs(Npro, label, distanceFile, Real, iFaceBalanced, proBalanced, stdBalanced,proList, NedgeRead, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames, plist, distMetrics);
  
  deltaFile<<label<<'\t'<<distMetrics[1]<<'\t'<<distMetrics[4]<<'\t'<<0<<'\t'<<0<<'\t'<<0<<'\t'<<0<<'\t'<<0<<'\t'<<0<<endl;
  /*Now do it again for the knocked down list
    replace proList, constrainKD, nPro, nConstrain. 
   */
  ofstream kdMatFile;
  ofstream balancedFile;
  kdMatFile.precision(12);
  balancedFile.precision(12);
  sprintf(fname, "KD_%s_RBratio_vsOriginal.txt", fappend.c_str());//label is a string.
  kdMatFile.open(fname);
  sprintf(fname, "KD_%s_BalancedCopies_vsOriginal.txt", fappend.c_str());//label is a string.
  balancedFile.open(fname);
  /*First print out all proteins, and their values in the original solution*/
  init_matrix_file(balancedFile, kdMatFile, proList, proBalanced);


  vector<double> RealKD;
  vector<double> iFaceBalancedKD;//Balanced solution--interfaces on a protein are not necessarily the same conc.
  vector<double> proBalancedKD; //Balanced solution after each protein is assigned an average over all interfaces
  vector<double> stdBalancedKD; //stand. deviation within a protein of interface copies.
  vector<double> distMetricsKD; //return the CSD and JSD for the balanced vs real copies.
  if(loopAll==false){
    
    distanceFile<<label2 << " Real to R'\tReal to R''\tR' to R''\n";
    NproWithKD=Npro-nProKD;
    solve_balanced_from_inputs(NproWithKD, label2, distanceFile, RealKD, iFaceBalancedKD, proBalancedKD, stdBalancedKD,proListKD, NedgeRead, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames, plist, distMetricsKD);
    /*Now here measure the distances between BalancedKD and Balanced. Calculate CSD and JSD between their distributions, normalize them first and leave out the 
      knocked-down protein.
     */
    /*Compare the Distributions before and after the KD*/
    compare_knockdown_with_original(deltaFile,balancedFile,  kdMatFile, label2, proListKD, proKDID, proList, RealKD, iFaceBalancedKD, proBalancedKD, iFaceBalanced, proBalanced);
    
  }else{
    //for each protein

    ofstream tmpfile;//("tmpfile.txt");
    
    for(int j=0;j<proList.size();j++){
      /*CLEAR VECTOR ARRAYS SO THEY ARE OVERWRITTEN.*/
      reset_arrays(knockDownNames,proKDID, proListKD, RealKD, iFaceBalancedKD, proBalancedKD, stdBalancedKD);
      /*Read in protein name proList[j] */
      tmpfile.open("tmpfile.txt");
      tmpfile <<proList[j].name;
      tmpfile.close();
      knockDownFile.open("tmpfile.txt");
      read_knockdowns(knockDownFile, knockDownNames, Npro, proList, proKDID, nProKD);
      knockDownFile.close();
      /*GENERATE NEW NETWORKS WITH THIS PROTEIN LIST, MINUS KNOCKDOWNS*/
      create_new_prolist(proList, proListKD,  proKDID,  nProKD);
      label2 = "KD";
      for(int i=0;i<nProKD;i++)
	label2=label2+":"+knockDownNames[i];
      
      distanceFile<<label2 << " Real to R'\tReal to R''\tR' to R''\n";
      NproWithKD=Npro-nProKD;
      solve_balanced_from_inputs(NproWithKD, label2, distanceFile, RealKD, iFaceBalancedKD, proBalancedKD, stdBalancedKD,proListKD,  NedgeRead, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames, plist, distMetricsKD);

      /*Compare the Distributions before and after the KD*/
      compare_knockdown_with_original(deltaFile,balancedFile,  kdMatFile, label2, proListKD, proKDID, proList, RealKD, iFaceBalancedKD, proBalancedKD, iFaceBalanced, proBalanced);
      
      
    }//done knocking down each protein.

  }
  cout<<" FINISHED ALL STOICHIOMETRIC BALANCE CALCULATIONS . Successful. "<<endl;
  distanceFile.close();
}


