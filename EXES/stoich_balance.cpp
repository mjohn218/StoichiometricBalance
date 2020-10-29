/*

./stoich_balance.exe parms.inp prolist_abund_constrain.inp edgelist.inp (CellType_AllProteinCopyNumber_Distribution.txt)  >std_out.out
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
struct MD_Timer looptime;
struct MD_Timer qptime;
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
    

  int Npro=plist.nwhole;

  string *genid=new string[Npro];
  double *abund=new double[Npro];
  vector<int> constrain;
  
  vector<ProteinClass> proList;
  proList.reserve(Npro);
  int nConstrain=0;
  string text;
  int t=0;
  
  /*READ IN LIST OF PROTEINS AND KNOWN CONCENTRATIONS */
  ifstream expfile(argv[2]);
  read_proteins_from_file( expfile,  plist,  genid,  abund,  proList,  nConstrain, constrain);

  
  /*READ IN NETWORK OF INTERACTIONS BETWEEN PROTEIN INTERFACES */
  ifstream iinfile(argv[3]);
  /*Read P1 I1 P2 I2*/
  vector<string> edge1Interface;
  vector<string> edge2Interface;
  vector<string> edge1Protein;
  vector<string> edge2Protein;
  vector<string> stoichiometry;
  
  read_edges_from_file(iinfile, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames);
  int NedgeRead=edge1Interface.size();
  
  /*REMOVE DUPLICATES FROM THE INTERFACE INTERACTION NETWORK, THESE CAN OCCUR IF YOU READ PAST THE END OF FILE.*/
  remove_duplicate_edges(NedgeRead, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry);

  /*IF CALCULATING P-VALUE, WILL READ IN A PROTEIN COPY_NUMBER DISTRIBUTION TO SAMPLE FROM*/
  vector<double> bin;
  vector<double> cdf;
  int nruns=plist.nruns_rand; //Number of random runs
  if(argc==5){  
    read_copynumber_dist(argv[4], bin,cdf, plist);
  }else{
    if(plist.flagread != 1)
      cout<<"For P-VALUE: Using Copy Numbers shuffled from the real COPY NUMBERS in Protein file: "<<argv[2]<<endl;
    else{
      nruns=0;
      cout<<" NO CELL TYPE DISTRIBUTION FILE READ IN. NOT PERFORMING STOICHIOMETRIC BALANCE FOR ANY RANDOMIZED COPY NUMBERS--NO P-VALUE CALCULATION."<<endl;
    }
  }
  
  /*
    POPULATE THE proLIST network, based on the proteins read in from file, matched to all the edges read in.
    SOME EDGES CONTAIN PROTEINS NOT IN THE INPUT FILE, THESE WILL BE SKIPPED.
    Construct the interface Network.
  */
  int nTotalInterfaces=0;
  int maxInterfacePartners=0;
  int nTotalEdges=0;
  
  construct_network_from_inputs( NedgeRead,  nTotalInterfaces,  maxInterfacePartners,  nTotalEdges, proList, edge1Interface, edge2Interface, edge1Protein, edge2Protein, stoichiometry, filenames);

  string *iNames=new string[nTotalInterfaces];
  int **Speclist=new int*[nTotalInterfaces];
  double **Stoichlist=new double*[nTotalInterfaces];
  for(i=0;i<nTotalInterfaces;i++){
      Speclist[i]=new int[maxInterfacePartners];
      Stoichlist[i]=new double[maxInterfacePartners];
  }
  int *numpartners=new int[nTotalInterfaces];
  Protein *wholep=new Protein[proList.size()];
  int Nif=nTotalInterfaces;
  int *p_home=new int[Nif];//this reverses and tells you what protein a given interface belongs to
  
  /*USE THE PROLIST NETWORK TO FORMAT THE INTERFACE ARRAYS FED INTO THE QP-SOLVERS.
   */
  construct_formatted_interface_arrays(nTotalInterfaces, proList,  wholep,  numpartners,  iNames, Speclist,  Stoichlist, p_home);
  
  int Ncomplex=0;//Need to account for extra edges in IIN
  string *edgenames=new string[nTotalEdges];
  define_edgenames( Ncomplex,  Nif,  numpartners,  Speclist,  edgenames,  nTotalEdges, iNames);
   
  
  /*BUILD THE A MATRIX, USED IN THE QP-SOLVERS. 
    A[i*Nif+j]: i is column, j is the row, each row is a different interface, each column a different complex.
  */
  
  double *A=new double[Nif*Ncomplex];
  for(i=0;i<Nif*Ncomplex;i++)
    A[i]=0;
    
  build_Amatrix_with_Stoich(Nif, A, numpartners, Speclist, Stoichlist);
  

  double *indivconc=new double[Nif];
  double *complexconc=new double[Ncomplex];
  double *complexmut=new double[Ncomplex];
  double *complexopt=new double[Ncomplex];
  double *curr_indiv=new double[Nif];


  
  /*
    INITIALIZE ALL ARRAYS FOR THE QP-SOLVE. THEY ARE DEFINED BASED ON THE NETWORK STRUCTURE, AND THE NUMBER OF CONSTRAINED COPY NUMBERS
    PARAMETER FOR SOLVING IS ALPHA: ENCODED HERE AS alphaValue
  */
  double alphaValue=plist.ascale;
  double nProCopies;
  
   /*Memory allocation for QP calls*/
  double *Q=new double[Ncomplex*Ncomplex];
  double *c=new double[Ncomplex];
  int nx=Ncomplex;
  /*bounds, should be zero*/
  double  *xupp=new double[nx];//0;//[] = { 20,   0 };
  char   *ixupp=new char[nx];//0;//[] = {  1,   0 };
  
  double  *xlow=new double[nx];//0;//[] = {  0,   0 };
  char   *ixlow=new char[nx];//0;//[] = {  1,   1 };
  const int nnzQ = Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2

  int    *irowQ=new int[nnzQ];// = {  0,   1,   1 }; 
  int    *jcolQ=new int[nnzQ];//[] = {  0,   0,   1 };
  double   * dQ=new double[nnzQ];//[] = {  8,   2,  10 };
  const int mz=Ncomplex;
  double *clow=new double[mz];
  char  *iclow=new char[mz];
  
  double *cupp=new double[mz];
  char  *icupp=new char[mz];
  const int nnzC = Ncomplex;
  int   *irowC=new int[nnzC];
  int   *jcolC=new int[nnzC];
  double   *dC=new double[nnzC];
  
  /*Construct same interfaces constraint matrix for this IIN*/
  double *H=new double[Nif*Nif];//this is size ninterfacexninterface
  /*construct matrix for weighting which protein concentrations are fixed*/
  double *ZA=new double[Nif*Ncomplex];
  double *Q2=new double[Nif*Ncomplex];
  
  initialize_qpMatrices( Ncomplex,  Nif,  nConstrain, plist, wholep, constrain, A, p_home,  xupp,  ixupp,  xlow,  ixlow,  irowQ,  jcolQ,   clow,  iclow,  cupp,  icupp,  irowC,  jcolC,  dC,  H,  ZA);
  
  /*****PART 2*****/
  /*Set up output files*/
  ofstream CSDfile;
  CSDfile.open("ChiSquaredDistances.out", ios::out | ios::trunc);
  ofstream ENTRfile;
  ENTRfile.open("JensenShannonDistances.out", ios::out | ios::trunc);
  
  CSDfile << "Run\tR to R'\tR to R''\tR' to R''\tReal to R (proteins)\n";
  ENTRfile << "Run\tR to R'\tR to R''\tR' to R''\tReal to R (proteins)\n";

  double *R = new double[Nif];
  double *Rp = new double[Nif];
  double *Rpp = new double[Nif];
  double *ivars = new double[Npro];
  double *imean = new double[Npro];
  double *CSD = new double[4];
  double *ENTR = new double[4];
  
  //Write variances of iface copy numbers for Rp on each protein to a file
  ofstream varfile;
  ofstream meanfile;
  varfile.open("Vars.out", ios::out | ios::trunc);
  meanfile.open("Means.out", ios::out | ios::trunc);
  varfile << "Run";
  for(i=0;i<Npro;i++){
    varfile << "\t" << genid[i];
    meanfile<<"\t"<<genid[i];
  }
  varfile<<endl;
  meanfile<<endl;
  //Also want to write R, R', and R'' for real copy numbers
  
  ofstream solvefile;
  solvefile.precision(12);
  solvefile.open("RealCopyNumbers_vsBalancedCopyNumbers.out", ios::out | ios::trunc);
  solvefile << "Protein\tInterface\tR\tR'\tR''\tStd for R'"<<endl;

  /***** Part 3 *****/
  initialize_timer(&looptime);
  initialize_timer(&qptime);
  start_timer(&looptime);

  /*Now that everything has been defined, do qp_solve*/
  /*We want to do this both with the real and random protein concentrations
    The first run (nc=0) will do the real, observed copy numbers for the proteins.
   */

  double pJSD=0;
  double pCSD=0;
  double realCSD,realJSD;
  
  int check=0;//Use this if some proteins not constrained
  double *initCopies = new double[Npro];
  //Set initCopiess to abund for first run
  for(i=0;i<Npro;i++)
    initCopies[i]=abund[i];
  cout<<"----------"<<endl;
  if(nruns>0)cout <<"NRUNS RANDOM COPY NUMBERS, to evaluate p-value: "<<nruns<<endl;

  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"---------------------PERFORM STOICHIOMETRIC BALANCE-------------------"<<endl;
  for(int nc=0;nc<nruns+1;nc++){
      start_timer(&qptime);
      
    if(nc>0){//Use random copy numbers
	cout <<"-------------------"<<endl;
	cout <<"Solving balance for the RANDOMIZED COPY NUMBERS,  to define a p-value. Run Number: "<<nc<<endl;
	if(plist.flagread==1)//Substituting flagread to determine what randomization method to use
	    sample_dist(Npro, bin,  cdf, initCopies);
	else
	    //Should only use shuffle when not constraining proteins!
	    shuffle_abund(Npro,abund,initCopies,nConstrain,constrain);
	/*cout <<"INITIAL COPY NUMBERS : "<<endl;
	//Set "abund" to initCopiess
	for(i=0;i<Npro;i++){
	    cout <<initCopies[i]<<endl;
	}//abund[i]=1.0*initCopies[i];
	*/
      //}
    }else{
	cout <<"-------------------"<<endl;
	cout <<"--------SOLVING BALANCE FOR THE REAL, OBSERVED PROTEIN COPY NUMBERS------"<<endl;
    }
    //cout << "Begin qp_solve: " << endl;

    //Make "R" from initCopies using p_home
    for(i=0;i<Nif;i++){
      R[i]=initCopies[p_home[i]];
    }

    //replacing "abund" with "initCopies" since initCopies is going to be variable
    
    qp_solve(Npro, Nif, Ncomplex,  wholep, indivconc, complexconc, A, plist, nConstrain, constrain, initCopies, Q, c, xlow, xupp, ixlow, ixupp, irowQ, jcolQ, dQ, clow, cupp, iclow, icupp, irowC, jcolC, dC, H, ZA, Q2, alphaValue,p_home);      
      
    stop_timer(&qptime);    
    //cout <<"QPsolve time: "<<timer_duration(qptime)<<endl;
      
    matmultiply(Nif,  Ncomplex, indivconc, complexconc, A);//Forward solve for Rp
       
    for(i=0;i<Nif;i++){
      Rp[i]=indivconc[i];
      //cout <<"i: "<<iNames[i]<<" "<<indivconc[i]<<endl;
    }
    if(nc==0){
	
	cout<<"-------------------------------"<<endl;
	cout <<" ------COPIES OF EACH COMPLEX FOR REAL COPY NUMBERS--------------"<<endl;
	for(i=0;i<Ncomplex;i++){
	    cout <<"edge: "<<i<<" "<<edgenames[i]<<" conc: "<<complexconc[i]<<endl;
	}
	cout<<"-------------------------------"<<endl;
	/*cout <<"--------Solved for interface copy numbers-----"<<endl;
	for(i=0;i<Nif;i++){
	    cout <<"i: "<<iNames[i]<<" "<<indivconc[i]<<endl;
	    }*/
    }
    //We can now get the variances using Rp
    get_variances(Npro, Nif, wholep, Rp, ivars, imean);
    //Now write to file
    if(nc==0){
      varfile<<"Real\t";
      meanfile<<"Real\t";
    }
    else{
      varfile<<"Random" << nc << "\t";
      meanfile<<"Random" << nc << "\t";
    }
    for(i=0;i<Npro;i++){
      varfile << ivars[i] << "\t";
      meanfile << imean[i] << "\t";
    }
    varfile << endl;
    meanfile << endl;


    even_protein(Npro, wholep, indivconc, nProCopies);//This changes indivconc by averaging concentrations for each protein. Therefore, I should save indivconc beforehand.
    for(i=0;i<Nif;i++){
      Rpp[i]=indivconc[i];
    }
      
    int i1, p1;
    
    /*First run is the Real Observed Copy numbers, write out full solution*/
    if(nc==0){
      for(p1=0;p1<proList.size();p1++){
	    for(i1=0;i1<proList[p1].InterfaceList.size();i1++){
		int index=proList[p1].InterfaceList[i1].globalIndex;
		solvefile<<proList[p1].name<<"\t"<<proList[p1].InterfaceList[i1].name<<"\t"<<R[index]<<'\t'<<Rp[index]<<'\t'<<Rpp[index]<<'\t'<<sqrt(ivars[p1])<<endl;
	    }//all interfaces on pro p1
	}
    }

    /*For looking at stats, want to look at constrained proteins only.*/
    if(nConstrain<Npro){
      for(j=0;j<Nif;j++){
	for(i=0;i<nConstrain;i++){
	  if(p_home[j]==constrain[i])
	    check=1;
	}
	if(check==0){
	  Rp[j]=R[j];
	  Rpp[j]=R[j]; //Making these equal means distance of zero
	}
	check=0;
      }

    }
    
    /*Next find distances from R to R', R to R'', and R' to R''.
Also want to look at distance between real copy numbers and random copy numbers*/

    CSD[0]=calc_chisquare(Nif,R,Rp);
    CSD[1]=calc_chisquare(Nif,R,Rpp);
    CSD[2]=calc_chisquare(Nif,Rp,Rpp);
    CSD[3]=calc_chisquare(Npro,abund,initCopies);

    //Calculating Jensen Shannon distance, rather than KL divergence
    ENTR[0]=calc_dist_jensenshannon(Nif,R,Rp,genid);
    ENTR[1]=calc_dist_jensenshannon(Nif,R,Rpp,genid);
    ENTR[2]=calc_dist_jensenshannon(Nif,Rp,Rpp,genid);
    ENTR[3]=calc_dist_jensenshannon(Npro,initCopies,abund,genid);
    
    //Write out entries to file
    if(nc==0){
      CSDfile << "Real\t" << CSD[0] << "\t" << CSD[1] << "\t" << CSD[2] << endl;
      ENTRfile << "Real\t" << ENTR[0] << "\t" << ENTR[1] << "\t" << ENTR[2] << endl;
    }
    else{
      CSDfile << "Random" << nc <<"\t" << CSD[0] << "\t" << CSD[1] << "\t" << CSD[2] << "\t" << CSD[3] << endl;
      ENTRfile << "Random" << nc <<"\t" << ENTR[0] << "\t" << ENTR[1] << "\t" << ENTR[2] << "\t" << ENTR[3] << endl;
    }

    //P-value calculations for R to R'
    if(nc==0){
      realCSD=CSD[0];
      realJSD=ENTR[0];
    }
    else{
      if(CSD[0]<realCSD)
	pCSD++;
      if(ENTR[0]<realJSD)
	pJSD++;
    }
      

  }//end looping over starting protein concentrations

  stop_timer(&looptime);    
  cout <<"Full time: "<<timer_duration(looptime)<<endl;

  varfile.close();
  meanfile.close();
  solvefile.close();

  pCSD/=(1.0*nruns);
  pJSD/=(1.0*nruns);

  CSDfile << "P-value for R to R': " << pCSD << endl;
  ENTRfile << "P-value for R to R': " << pJSD << endl;
  

  CSDfile.close();
  ENTRfile.close();
  
}


