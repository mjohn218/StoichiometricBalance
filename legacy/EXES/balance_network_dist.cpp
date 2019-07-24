/*

./balance_network_dist.exe parms.inp node_abund.inp edgelist.inp ref_IIN.inp Probs_abundance.txt [constraints.inp]

Oct 2017 (final version)

Inputs
- parms.inp                set of parameters
- node_abund.inp           abundance of proteins
- edgelist.inp             list of protein-protein interactions
- ref_IIN.inp              The IIN (interface interaction network)
- Probs_abundance.txt      Protein concentration distribution to be used to sample random copy numbers
- constraints.inp          Optional input selecting which proteins have constrained concentrations. Default is that all proteins are constrained with numbers from node_abund.inp

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

#include "md_timer.h"
#include "pro_classes.h"
#include "constrainParms.h"
#include "metrics_class.h"
#include "matmultiply.h"
#include "qp_solvers.h"

#include "read_proinput.h"
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

void sample_dist(int npro, int *bin, double *cmf, double *cdf, double *rval);
void shuffle_abund(int npro, double *abund, double *rval, int Nconstrain, int*constrain);
void get_variances(int Npro, int Nif, Protein *wholep, double *R, double *ivars, double *imean);

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
  //cout <<"seed: "<<seed<<endl;
  
  srand(seed);
  
  ifstream parmfile(argv[1]);
  constrainParms plist;
  read_balance_parms(parmfile, plist);
  write_parms(plist);

    
  /*read in the list of proteins in the network and their known concentration*/
  int Npro=plist.nwhole;
  ppidata *ppi=new ppidata[Npro];
  string *genid=new string[Npro];
  double *abund=new double[Npro];
  
  
  ifstream expfile(argv[2]);
  for(i=0;i<Npro;i++){
    expfile>>genid[i]>>abund[i];
    //cout <<genid[i]<<' '<<abund[i]<<endl;
  }
  cout <<"read in abundances "<<endl;
  /*read in the network of protein protein interactions*/
  int Nedge=plist.Nedge;
  ifstream edgefile(argv[3]);
  string *e1=new string[Nedge];
  string *e2=new string[Nedge];
  string ig;
  edgefile.ignore(400,'\n');//skip first line
  for(i=0;i<Nedge;i++){
    edgefile >>e1[i]>>ig>>e2[i];
    //cout <<e1[i]<<' '<<e2[i]<<endl;
  }
  cout <<"read in edges "<<endl;
  double *Padj=new double[Npro*Npro];
  for(i=0;i<Npro*Npro;i++)
    Padj[i]=0;

  for(i=0;i<Npro;i++)
    ppi[i].nppartner=0;
  
  int p;
  int p1, p2;
  int *e1num=new int[Nedge];
  int *e2num=new int[Nedge];
  int t=0;

  for(i=0;i<Nedge;i++){
    /*its possible there is an edge in the network with a node not included */
    p1=-1;
    p2=-1;
    for(j=0;j<Npro;j++){
      if(e1[i]==genid[j])
	p1=j;
      if(e2[i]==genid[j])
	p2=j;
    }
    if(p1==-1 || p2==-1){
      cout <<"unmatched edge: "<<i<<" skipping..." <<endl;

    }
    else{
      e1num[t]=p1;
      e2num[t]=p2;
      t++;
      /*Protein-Protein adjacency matrix*/
      Padj[p1*Npro+p2]=1;
      if(p1==p2){
	//self interaction
	//cout <<"self interaction: "<<p1<<' '<<genid[p1]<<endl;
	Padj[p2*Npro+p1]=2;
	p=ppi[p1].nppartner;
	ppi[p1].pplist[p]=p1;
	ppi[p1].nppartner++;
      }
      else{
	p=ppi[p1].nppartner;
	ppi[p1].pplist[p]=p2;
	ppi[p1].nppartner++;
	p=ppi[p2].nppartner;
	ppi[p2].pplist[p]=p1;
	ppi[p2].nppartner++;
      }
    }
  }
  Nedge=t;
  //cout <<"total protein interactions: "<<Nedge<<endl;
  //  for(i=0;i<Npro;i++)
    //   cout <<i<<' '<<genid[i]<<' '<<ppi[i].nppartner<<endl;
  //cout <<"set Protein-Protein adjacency matrix "<<endl;
  /*read in data on proteins levels to constrain*/
  int Nconstrain;
  if(argc==7)
    Nconstrain=plist.Nconstrain;
  else
    Nconstrain=Npro;

  int *constrain=new int[Nconstrain];
  
  /*assign a unique interface to each protein*/
  int Nif=Nedge*2.0;//this is actually the maximum, it will be less if there are self interactions
  Protein *wholep=new Protein[Npro];
  int *p_home=new int[Nif];//this reverses and tells you what protein a given interface belongs to
  t=0;
  int nint;
  /*Be careful with self interactions, for now, just also make them have a distinct interface, but it
   binds to itself, otherwise you create a chain, head to tail , head to tail*/

  /*READ IN THE ACTUAL INTERFACE DISTRIBUTION*/
  int **Speclist=new int*[Nedge*2];
  for(i=0;i<Nedge*2;i++)
    Speclist[i]=new int[MAXP];
  int *numpartners=new int[Nedge*2];//this is just maximum number of interfaces

  //if(plist.flagread==1){
  ifstream reffile(argv[4]);
  cout <<"reading in IIN from file "<<argv[4]<<endl;
  read_ref_ProIIN(Npro, wholep, plist.ninterface, reffile, numpartners, Speclist, p_home);
  cout <<"number of interfaces: "<<plist.ninterface<<endl;
  Nif=plist.ninterface;

  ifstream distfile(argv[5]);
  /*Read in conc distribution*/
  int MAXBINS=5000;
  int *bin = new int[MAXBINS];
  double *cmf = new double[MAXBINS];
  double *cdf = new double[MAXBINS];

  t=1;
  bin[0]=0;
  cmf[0]=0;
  cdf[0]=0;
  while(distfile >> bin[t] >> cmf[t] >> cdf[t])
    t++;

  distfile.close();
     
  /*Create the matrix A*/
  int Ncomplex;
  
  Ncomplex=0; //Need to account for extra edges in IIN
  for(i=0;i<Nif;i++){
    for(j=0;j<numpartners[i];j++){
      if(Speclist[i][j]>=i){//Avoid counting an edge twice
	Ncomplex++;
      }
    }
  }
   
  cout << "IIN edges: " << Ncomplex << endl;
  double *A=new double[Nif*Ncomplex];
  // cout <<"matrix size: "<<Nif*Ncomplex<<endl;
  for(i=0;i<Nif*Ncomplex;i++)
    A[i]=0;
    
  build_Amatrix(Nif, A, numpartners, Speclist);
  
  /*A[i*Nif+j] i is column, j is the row, each row is a different interface, each column a different complex.*/
    
   

  /*Establish constraints on the number of interfaces, based on protein home */
  double *indivconc=new double[Nif];
  double *complexconc=new double[Ncomplex];
  double *complexmut=new double[Ncomplex];
  double *complexopt=new double[Ncomplex];
  double *curr_indiv=new double[Nif];
  
  /*The total number of protein complexes is half the total number of interfaces when perfectly matched*/
  /*guess the number of protein complexes, total ifaces is numiface. Total complex numbers can change however when protein
   concentrations change*/

  int numstart;//=conc_start*V;
  int ni;

    
  double ascale=plist.ascale;
  
  double *randvals=new double[Nif];
  double delfit, rnum;
  char fname[200];
  
  double npfinal;
  char onec[200];
  char longc[400];
  double *proconc=new double[Npro];
  double entropy, entropy_rand;
  double msd, msd_rand;
  
    
  /*Memory allocation for QP calls*/
  double *Q=new double[Ncomplex*Ncomplex];
  double *c=new double[Ncomplex];
  int nx=Ncomplex;
  /*bounds, should be zero*/
  double  *xupp=new double[nx];//0;//[] = { 20,   0 };
  char   *ixupp=new char[nx];//0;//[] = {  1,   0 };
  
  double  *xlow=new double[nx];//0;//[] = {  0,   0 };
  char   *ixlow=new char[nx];//0;//[] = {  1,   1 };
  for(i=0;i<nx;i++){
    xupp[i]=0;
    ixupp[i]=0;
    xlow[i]=0;
    ixlow[i]=0;
  }
  const int nnzQ = Ncomplex*(Ncomplex+1)/2;//Ncomplex choose 2

  int    *irowQ=new int[nnzQ];// = {  0,   1,   1 }; 
  int    *jcolQ=new int[nnzQ];//[] = {  0,   0,   1 };
  double   * dQ=new double[nnzQ];//[] = {  8,   2,  10 };
  t=0;
  for(i=0;i<Ncomplex;i++){
    for(j=i;j<Ncomplex;j++){
      jcolQ[t]=i;
      irowQ[t]=j;
      //dQ[t]=Q[i*Ncomplex+j];
      t++;
    }
  }
  const int mz=Ncomplex;
  double *clow=new double[mz];
  char  *iclow=new char[mz];
  
  double *cupp=new double[mz];
  char  *icupp=new char[mz];
  for(i=0;i<mz;i++){
    clow[i]=plist.min_complex;
    iclow[i]=1;
    cupp[i]=0;
    icupp[i]=0;
  }
  const int nnzC = Ncomplex;
  int   *irowC=new int[nnzC];
  int   *jcolC=new int[nnzC];
  double   *dC=new double[nnzC];
  for(i=0;i<Ncomplex;i++){
    irowC[i]=i;
    jcolC[i]=i;
    dC[i]=1.0;
  }


  /*Construct same interfaces constraint matrix for this IIN*/
  double *H=new double[Nif*Nif];//this is size ninterfacexninterface
  /*construct matrix for weighting which protein concentrations are fixed*/
  double *Z=new double[Nif*Nif];
  double *ZA=new double[Nif*Ncomplex];
  double *Q2=new double[Nif*Ncomplex];
  double cof;
  for(i=0;i<Nif;i++){
    p1=p_home[i];
    nint=wholep[p1].ninterface;
    cof=1.0;
    H[i*Nif+i]=(nint-1)/cof;
    //cout <<"Hi: "<<H[i*Nif+i]<<" nint: "<<nint<<endl;
    
    for(j=i+1;j<Nif;j++){
      p2=p_home[j];
      if(p2==p1){
	H[i*Nif+j]=-1.0/cof;
	H[j*Nif+i]=-1.0/cof;
      }
      else{
	H[i*Nif+j]=0;
	H[j*Nif+i]=0;
	
      }
    }
  }


  /*Begin looping over distinct sets of constraints*/
  initialize_timer(&looptime);
  initialize_timer(&qptime);
  start_timer(&looptime);
  Metrics mets;
  Metrics rndmets;
  double rnd_rmsd_avg;
  double avg_degree=2.0*Ncomplex/(1.0*Nif);
  double constrn_fit;
  //cout <<"Nruns random: "<<plist.nruns_rand<<endl;
  
  //Start out Z as zero matrix
  for(i=0;i<Nif;i++){      
    Z[i*Nif+i]=0;
  }

  int ntmp;
  if(argc==7){
    ifstream constrainfile(argv[6]);
    for(i=0;i<Nconstrain;i++){
      constrainfile >> constrain[i];
      for(j=0;j<wholep[constrain[i]].ninterface;j++){
	ntmp=wholep[constrain[i]].valiface[j];
	Z[ntmp*Nif+ntmp]=1;
      }
    }
  }
  else{
    for(i=0;i<Npro;i++)
      constrain[i]=i;
    //Set Z to identity matrix
    for(i=0;i<Nif;i++){      
      Z[i*Nif+i]=1;
    }
  }
  
  

  /*define Z*A matrix*/
  for(i=0;i<Nif;i++){
    if(Z[i*Nif+i]==1){
      //copy A into ZA
      for(j=0;j<Ncomplex;j++){
	ZA[j*Nif+i]=A[j*Nif+i];
      }
    }else{
      for(j=0;j<Ncomplex;j++)
	ZA[j*Nif+i]=0;
    }
  }

  if(plist.flagread==1)
    cout<<"Using randomly sampled CNs\n";
  else
    cout<<"Using CNs shuffled from the real CNs\n";
  
  /*****PART 2*****/

  /*Set up output files*/

  ofstream CSDfile;
  CSDfile.open("CSD.out", ios::out | ios::trunc);
  ofstream ENTRfile;
  ENTRfile.open("JSD.out", ios::out | ios::trunc);
  
  CSDfile << "Run\tR to R'\tR to R''\tR' to R''\tReal to R (proteins)\n";
  ENTRfile << "Run\tR to R'\tR to R''\tR' to R''\tReal to R (proteins)\n";

  double *R = new double[Nif];
  double *Rp = new double[Nif];
  double *Rpp = new double[Nif];
  double *ivars = new double[Npro];
  double *imean = new double[Npro];
  
  int nruns=plist.nruns_rand; //Number of random runs
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
  solvefile.open("SolvedProteins.out", ios::out | ios::trunc);
  solvefile << "Protein\tA\tA'\tA''\tvar for A'"<<endl;

  /***** Part 3 *****/
  
  /*Now that everything has been defined, do qp_solve*/
  /*We want to do this both with the real and random protein concentrations*/

  double pJSD=0;
  double pCSD=0;
  double realCSD,realJSD;
  
  int check=0;//Use this if some proteins not constrained
  double *rval = new double[Npro];
  //Set rvals to abund for first run
  for(i=0;i<Npro;i++)
    rval[i]=abund[i];

  for(int nc=0;nc<nruns+1;nc++){
    start_timer(&qptime);
    if(nc>0){//Use random copy numbers
      if(plist.flagread==1)//Substituting flagread to determine what randomization method to use
	sample_dist(Npro, bin, cmf, cdf, rval);
      else
	//Should only use shuffle when not constraining proteins!
	shuffle_abund(Npro,abund,rval,Nconstrain,constrain);
      //Set "abund" to rvals
      //for(i=0;i<Npro;i++){
      //abund[i]=1.0*rval[i];
      //}
    }
    //cout << "Begin qp_solve: " << endl;

    //Make "R" from rval using p_home
    for(i=0;i<Nif;i++){
      R[i]=rval[p_home[i]];
    }

    //replacing "abund" with "rval" since rval is going to be variable
    
    qp_solve(Npro, Nif, Ncomplex, numstart, wholep, indivconc, complexconc, A, plist, Nconstrain, constrain, rval, Q, c, xlow, xupp, ixlow, ixupp, irowQ, jcolQ, dQ, clow, cupp, iclow, icupp, irowC, jcolC, dC, H, ZA, Q2, ascale,p_home);      
      
    stop_timer(&qptime);    
    //cout <<"QPsolve time: "<<timer_duration(qptime)<<endl;
      
    matmultiply(Nif,  Ncomplex, indivconc, complexconc, A);//Forward solve for Rp
      
    for(i=0;i<Nif;i++){
      Rp[i]=indivconc[i];
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


    even_protein(Npro, wholep, indivconc, npfinal);//This changes indivconc by averaging concentrations for each protein. Therefore, I should save indivconc beforehand.
    for(i=0;i<Nif;i++){
      Rpp[i]=indivconc[i];
    }
      
    int i1;
      //      cout <<"Protein concentrations, averaged "<<endl;
    for(i=0;i<Npro;i++){
      i1=wholep[i].valiface[0];
      proconc[i]=indivconc[i1];
    
      //cout <<proconc[i]<<endl;
    }

    //If real, write out the solutions
    if(nc==0){
      for(i=0;i<Nif;i++){
	solvefile << p_home[i] << "\t" << R[i] << "\t"<< Rp[i] << "\t" << Rpp[i] << "\t" << ivars[p_home[i]] << endl;
      }

    }

    /*For looking at stats, want to look at constrained proteins only.*/
    if(Nconstrain<Npro){
      for(j=0;j<Nif;j++){
	for(i=0;i<Nconstrain;i++){
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
    CSD[3]=calc_chisquare(Npro,abund,rval);

    //Calculating Jensen Shannon distance, rather than KL divergence
    ENTR[0]=calc_dist_jensenshannon(Nif,R,Rp,genid);
    ENTR[1]=calc_dist_jensenshannon(Nif,R,Rpp,genid);
    ENTR[2]=calc_dist_jensenshannon(Nif,Rp,Rpp,genid);
    ENTR[3]=calc_dist_jensenshannon(Npro,rval,abund,genid);
    
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

void sample_dist(int npro, int *bin, double *cmf, double *cdf, double *rval)
{
  /*This function samples protein distributions from an input file Probs_abundance.out*/

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  timeval tim;
  gettimeofday(&tim,0);
  double t1=tim.tv_sec + tim.tv_usec;
  int seed=t1;
  gsl_rng_set(r, t1);
  
  //int npro = atoi(argv[1]);
  
  int i;
  //int *rval = new int[npro];
  double rnum;
  int t;
  for(i=0;i<npro;i++){

    t=0;
    rnum=gsl_rng_uniform(r);
    //cout << "rnum is: " << rnum << endl;
    while(rnum>cdf[t])
      t++;

    t--;
    double fract=(cdf[t+1]-rnum)/(cdf[t+1]-cdf[t]);
    rval[i]=int(fract*(bin[t+1]-bin[t])+bin[t]);
  }
  //infile.close();
  
  //ofstream outfile;
  //outfile.open("protnums_rand.inp",ios::out|ios::trunc);
  //for(i=0;i<npro;i++)
  //outfile << rval[i] << endl;

  //outfile.close();

  gsl_rng_free(r);
  
}
void shuffle_abund(int npro, double *abund, double *rval, int Nconstrain, int *constrain)
{
  //Shuffle the copy number values for each protein.

  //Set random generator seed
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc(T);

  timeval tim;
  gettimeofday(&tim,0);
  double t1=tim.tv_sec + tim.tv_usec;
  int seed=t1;
  gsl_rng_set(r, t1);

  //First set rval to abund
  int i,j,i2,j2;
  double temp;
  for(i=0;i<npro;i++){
    rval[i]=abund[i];
  }
  
  double rnum;
  /*for(i = npro-1; i>0; i--){
    rnum=npro*gsl_rng_uniform(r);
    j=int(rnum);
    //Swap a[i] and a[j]
    temp = rval[i];
    rval[i]=rval[j];
    rval[j]=temp;

    }*/

  //Only want to shuffle proteins that we are constraining
  for(i = Nconstrain-1; i>0; i--){
    rnum=Nconstrain*gsl_rng_uniform(r);
    j=int(rnum);
    //Swap a[i2] and a[j2];
    j2=constrain[j];
    i2=constrain[i];
    temp = rval[i2];
    rval[i2]=rval[j2];
    rval[j2]=temp;

  }

  gsl_rng_free(r);

}


void get_variances(int Npro, int Nif, Protein *wholep, double *R, double *ivars, double *imean)
{

  double sum1,sum2;
  int i,j;
  int i1,ni;

  //Calculate mean and variance
  
  
  for(i=0;i<Npro;i++){
    sum1=0;
    sum2=0;
    ni=wholep[i].ninterface;
    for(j=0;j<ni;j++){
      i1=wholep[i].valiface[j];
      cout << i1 << ": " << R[i1] << endl;
      sum1+=R[i1];
      sum2+=R[i1]*R[i1];
    }
    //Now divide
    sum1/=(1.0*ni); //This is the mean
    sum2/=(1.0*ni); //This is E[x^2]
    ivars[i]=sum2 - sum1*sum1; //This is the variance
    imean[i]=sum1;
  }

}

