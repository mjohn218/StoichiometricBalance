/*

./stoich_balance.exe parms.inp prolist_abund_constrain.inp edgelist.inp CellType_AllProteinCopyNumber_Distribution.txt  >std_out.out



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

class Interface{

public:
    string name;
    //string proteinHome;
    int globalIndex;
    int nPartners;
    //    vector<string> partnerList;
    vector<string> proPartnerName;
    vector<string> interfacePartnerName;
    vector<int> proPartnerIndex;
    vector<int> interfacePartnerIndex;//global index
    Interface(int index);//constructor
    vector<double> stoichiometry;
};
Interface::Interface(int index)
{
    globalIndex=index;
}
class ProteinClass{
public:
    int nInterfaces;
    string name;
    //    vector<string> iFaceNames;
    vector<Interface> InterfaceList;//has properties of the above class
    int nProPartners;
    //vector<string> proPartnerNames;
    int copies;
    bool isConstrained;
    int index;
    ProteinClass(int n_copies);//constructor
};
ProteinClass::ProteinClass(int n_copies)
{
    copies=n_copies;
}
 
void sample_dist(int npro, vector<double> bin, vector<double> cdf, double *rval);
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
  plist.max_complex=0;//initialize in case it is not read in.
  
  read_short_parms(parmfile, plist);
  //write_parms(plist);
  cout <<"------------------ Parameters -----------"<<endl;
  cout <<" N Proteins : "<<plist.nwhole<<endl;
  cout <<" Minimum number of Complexes assignable (>=0): "<<plist.min_complex<<endl;
  cout <<" Alpha value assigned: "<<plist.ascale <<". A value of 1 works well. Very small forces all interfaces on same protein to be equal copy numbers. "<<endl;
  cout <<" Number of times to repeat with randomized copy numbers, solely to evaluate p-value : "<<plist.nruns_rand<<endl;
  cout <<" Only if above is >0: initialize random copies from the Cell Type distribution, or by shuffling the observed copy numbers: "<<plist.flagread<<". (Set to 1 to sample from Cell Type distribution.) "<<endl;
  cout <<" Max Complexes: "<<plist.max_complex<<". Default=0 will allow max to infinity (unconstrained max value). "<<endl;
  //cout <<"WARNING: If Max Complexes is a finite value, one cannot sample randomized copy numbers for p-value calculations that exceed this, so it should accomodate the maximum protein copy numbers in the CellType_Distribution_CopyNumbers file. "<<endl;
    
  /*read in the list of proteins in the network and their known concentration*/
  int Npro=plist.nwhole;

  string *genid=new string[Npro];
  double *abund=new double[Npro];
  string *applyConstraint=new string[Npro];
  vector<ProteinClass> proList;
  
  ifstream expfile(argv[2]);
  int Nconstrain=0;
  string text;
  proList.reserve(Npro);
  cout<<"----------"<<endl;
  //  vector<int> constrain;//=new int[Nconstrain];
  expfile.ignore(400,'\n');//>>text>>text>>text;
  ProteinClass temp(-1);
  Interface itemp(0);
  int cnt=0;
  for(i=0;i<Npro;i++){
      /*
	Read in name of the protein, its copy numbers (use dummy value if unknown) and
	either Yes or No for applyConstraint to that protein. 
      */
      expfile>>genid[i]>>abund[i]>>applyConstraint[i];
      //cout <<"genid: "<<genid[i]<<" abund: "<<abund[i]<<" isConstrained?: "<<applyConstraint[i]<<endl;
      proList.push_back(temp);
      proList[i].name=genid[i];
      proList[i].copies=abund[i];
      proList[i].nInterfaces=0;//haven't added these yet, initialize to zero
      if(applyConstraint[i]=="Y"){
	  proList[i].isConstrained=true;
	  //constrain.push_back(i);
	  Nconstrain++;
      }else
	  proList[i].isConstrained=false;
      //add this protein to the List.
      proList[i].index=i;
      cnt++;
      //cout <<genid[i]<<' '<<abund[i]<<endl;
  }
  /*cout <<" -------------------- "<<endl;
  cout <<" PROTEINS COUNTED: "<<cnt<<" FROM INPUT FILE: "<<argv[2]<<'\t';
  
  if(cnt!=Npro){
      cout <<" ************************************ "<<endl;
      cout <<"********WARNING: READ IN DIFFERENT NUMBER OF PROTEINS " <<cnt<<" THAN LISTED IN : "<<argv[1]<<", which is: "<<Npro<<" *********"<<endl;
      cout <<" ************************************ "<<endl;
      }*/
  int t=0;
  int *constrain=new int[Nconstrain];
  for(i=0;i<Npro;i++){
      if(applyConstraint[i]=="Y"){
	  constrain[t]=i;
	  t++;
      }
  }
  
  //  cout <<"read in abundances "<<endl;
  cout <<" -------------------- "<<endl;
  cout <<" NUMBER OF CONSTRAINED PROTEIN ABUNDANCES : "<<Nconstrain<<". Calculated from file: "<<argv[2]<<endl;
  
  /*read in the network of protein protein interactions*/
  ifstream iinfile(argv[3]);
  /*Read P1 I1 P2 I2*/
  vector<string> edge1Interface;
  vector<string> edge2Interface;
  vector<string> edge1Protein;
  vector<string> edge2Protein;
  vector<string> stoichiometry;
  string pro1, pro2;
  string iface1, iface2, stoich;
  int NedgeRead=0;
  iinfile.ignore(600,'\n');//first line is titles of columns
  //  iinfile >>text>>text>>text>>text;
  while(!iinfile.eof()){
      iinfile >>pro1 >>iface1>>pro2>>iface2>>stoich;
      edge1Interface.push_back(iface1);
      edge2Interface.push_back(iface2);
      edge1Protein.push_back(pro1);
      edge2Protein.push_back(pro2);
      stoichiometry.push_back(stoich);
      NedgeRead++;
  }
  cout <<" -------------------- "<<endl;
  cout <<"NUMBER OF EDGES : "<<NedgeRead<<" read from file: "<<argv[3]<<endl;
  /*REMOVE DUPLICATES FROM THE INTERFACE INTERACTION NETWORK, THESE CAN OCCUR IF YOU READ PAST THE END OF FILE.*/
  bool duplicate=false;
  for(i=1;i<NedgeRead;i++){
      //look at all previous edges for duplicates
      for(int prev=0;prev<i;prev++){
	  duplicate=false;
	  if(edge1Protein[prev]==edge1Protein[i] && edge2Protein[prev]==edge2Protein[i] && edge1Interface[prev]==edge1Interface[i] && edge2Interface[prev]==edge2Interface[i])
	      {
		  cout <<"READ IN THE SAME EXACT INTERACTION TWICE, DELETING ONE OF THEM. Lines: "<<prev+2<<" and "<<i+2<<endl;
		  /*Replace the previous edge with the last edge in the list, then decrease the size of the array.*/
		  int last=NedgeRead-1;
		  edge1Protein[prev]=edge1Protein[last];
		  edge2Protein[prev]=edge2Protein[last];
		  edge1Interface[prev]=edge1Interface[last];
		  edge2Interface[prev]=edge2Interface[last];
		  stoichiometry[prev]=stoichiometry[last];
		  NedgeRead--;
		  duplicate=true;

	      }
	  if(duplicate==false){
	      //check the reversed A B interaction, but not if B A is same for this pair.
	      if(edge1Protein[prev]==edge2Protein[i] && edge2Protein[prev]==edge1Protein[i] && edge1Interface[prev]==edge2Interface[i] && edge2Interface[prev]==edge1Interface[i])
		  {
		      cout <<"READ IN THE SAME EXACT INTERACTION TWICE, (REVERSED A AND B COLUMNS). DELETING ONE OF THEM. Lines: "<<prev+2<<" and " <<i+2<<endl;
		      /*Replace the previous edge with the last edge in the list, then decrease the size of the array.*/
		      int last=NedgeRead-1;
		      edge1Protein[prev]=edge1Protein[last];
		      edge2Protein[prev]=edge2Protein[last];
		      edge1Interface[prev]=edge1Interface[last];
		      edge2Interface[prev]=edge2Interface[last];
		      stoichiometry[prev]=stoichiometry[last];
		      NedgeRead--;
		      
		  }
	  }
      }
  }
  cout <<" -------------------- "<<endl;
  cout <<"AFTER CHECKING FOR DUPLICATES, EDGES READ IN IS NOW: "<<NedgeRead<<endl;
  cout <<" -------------------- "<<endl;
  int p;
  int p1, p2;
  /*Match Proteins in the interface network against the proteins in the
    main protein file.
    Construct the interface Network.
  */
  int nTotalInterfaces=0;
  int maxInterfacePartners=0;
  int nTotalEdges=0;
  cout <<" -------------------- "<<endl;
  cout <<" MATCH EDGES BASED ON PROTEIN LIST FROM FILE "<<argv[2]<<endl;
  double *stoichvalue1=new double[NedgeRead];
  double *stoichvalue2=new double[NedgeRead];
  string delimiter=":";
  for(i=0;i<NedgeRead;i++){
    /*its possible there is an edge in the network with a node not included */
    p1=-1;
    p2=-1;
    for(j=0;j<proList.size();j++){
      if(edge1Protein[i]==proList[j].name)
	p1=j;
      if(edge2Protein[i]==proList[j].name)
	p2=j;
    }
    if(p1==-1 || p2==-1){
	cout <<"CANNOT FIND BOTH PROTEINS IN EDGE: "<<i<<" skipping... " <<edge1Protein[i]<<' '<<edge2Protein[i]<<endl;
    }else{
	nTotalEdges++;
	//cout <<"Current edge: "<<i<<" proteins: "<<p1<<' '<<p2<<' '<<endl;
	bool doesExist1=false;//does the interface exist yet?
	bool doesExist2=false;//does the interface exist yet?
	/*Use this protein interaction to construct the interface interaction network*/
	int currIndex;
	int localIndex;
	//cout <<"P1: "<<p1<<" Current number of interfaces: "<<proList[p1].InterfaceList.size()<<" Interface to match: "<<edge1Interface[i]<<endl;
	for(int m=0;m<proList[p1].InterfaceList.size();m++){
	    //does this interface exist? Otherwise add it.
	    if(proList[p1].InterfaceList[m].name==edge1Interface[i]){
		//this interface already exists
		doesExist1=true;
		localIndex=m;
		currIndex=proList[p1].InterfaceList[m].globalIndex;//which interface this is of all.
		
	    }
	}//looping over existing nInterfaces
	/*split the stoichiometry ratio into two values, stoichvalue1 and stoichvalue2*/
	
	string str=stoichiometry[i].substr(0, stoichiometry[i].find(delimiter));
	
	int start=stoichiometry[i].find(delimiter)+delimiter.length();
	stoichvalue1[i]=stod(str);
	string str2=stoichiometry[i].substr(start, stoichiometry[i].length());
	stoichvalue2[i]=stod(str2);
	//cout <<" first ratio value: "<<str<<" second value "<<str2<<" doubles: "<<stoichvalue1[i]<<' '<<stoichvalue2[i]<<endl;
	if(doesExist1==false){
	    
	    proList[p1].nInterfaces++;
	    localIndex=proList[p1].InterfaceList.size();//new interface on the protein, size before increasing
	    //create a new interface for this protein.
	    proList[p1].InterfaceList.push_back(itemp);//add this interface to the list for this protein.
	    proList[p1].InterfaceList[localIndex].name=edge1Interface[i];
	    //proList[p1].InterfaceList[localIndex].proteinHome=edge1Protein[i];
	    proList[p1].InterfaceList[localIndex].globalIndex=nTotalInterfaces;//index in full interface list
	    currIndex=nTotalInterfaces;
	   
	    nTotalInterfaces++;
	    //proList[p1].iFaceNames.push_back(edge1Interface[i]);//add this interface to the protein.
	    
	    
	}
	/**/
	proList[p1].InterfaceList[localIndex].proPartnerName.push_back(edge2Protein[i]);//add protein
	proList[p1].InterfaceList[localIndex].proPartnerIndex.push_back(p2);//add protein
	proList[p1].InterfaceList[localIndex].interfacePartnerName.push_back(edge2Interface[i]);//add interface
	proList[p1].InterfaceList[localIndex].stoichiometry.push_back(stoichvalue1[i]);
	int currsize=proList[p1].InterfaceList[localIndex].proPartnerName.size();
	if(currsize>maxInterfacePartners){
	    maxInterfacePartners=currsize;
	}

	/*Now do P2, add this interface and the interaction to the proteins*/
	//	cout <<"P2: "<<p2<<" Current number of interfaces: "<<proList[p2].InterfaceList.size()<<" Interface to match: "<<edge2Interface[i]<<endl;
	if(p1==p2 && edge1Interface[i]==edge2Interface[i]){
	    //this is a self interaction, don't add it again. 
	    //cout <<"SELF BINDING INTERFACE: "<<edge1Interface[i]<<' '<<edge2Interface[i]<<" on Protein: "<<p1<<endl;
	}else{
	    for(int m=0;m<proList[p2].InterfaceList.size();m++){
		//does this interface exist? Otherwise add it.
		if(proList[p2].InterfaceList[m].name==edge2Interface[i]){
		    //this interface already exists
		    doesExist2=true;
		    localIndex=m;
		    
		    currIndex=proList[p2].InterfaceList[m].globalIndex;//which interface this is of all.
		}
	    }//looping over existing nInterfaces
	    
	    
	    if(doesExist2==false){
		
		proList[p2].nInterfaces++;
		localIndex=proList[p2].InterfaceList.size();//new interface on the protein, current index
		//create a new interface for this protein.
		proList[p2].InterfaceList.push_back(itemp);//add this interface to the list for this protein.
		
		proList[p2].InterfaceList[localIndex].name=edge2Interface[i];
		//proList[p2].InterfaceList[localIndex].proteinHome=edge1Protein[i];
		proList[p2].InterfaceList[localIndex].globalIndex=nTotalInterfaces;//index in full interface list
		currIndex=nTotalInterfaces;
		nTotalInterfaces++;
		
	    }
	    /**/
	    proList[p2].InterfaceList[localIndex].proPartnerName.push_back(edge1Protein[i]);//add protein
	    proList[p2].InterfaceList[localIndex].proPartnerIndex.push_back(p1);//add protein
	    proList[p2].InterfaceList[localIndex].interfacePartnerName.push_back(edge1Interface[i]);//add interface
	    proList[p2].InterfaceList[localIndex].stoichiometry.push_back(stoichvalue2[i]);//add interface
	    currsize=proList[p2].InterfaceList[localIndex].proPartnerName.size();
	    
	    if(currsize>maxInterfacePartners){
		maxInterfacePartners=currsize;
	    }
	}//non self, add P2 interface
	
    }
  }//done looping over all edges in the input file.
  
  int i1;
  int nint;
  /*Now figure out globalIndex of each interface's partner interface*/
  for(p1=0;p1<proList.size();p1++){
	for(j=0;j<proList[p1].InterfaceList.size();j++){
	    i1=proList[p1].InterfaceList[j].globalIndex;
	    for(int k=0;k<proList[p1].InterfaceList[j].interfacePartnerName.size();k++){
		string i2name=proList[p1].InterfaceList[j].interfacePartnerName[k];
		int p2=proList[p1].InterfaceList[j].proPartnerIndex[k];
		//find the global index of this interface i2name on protein p2.
		int i2=-1;
		for(int s=0;s<proList[p2].InterfaceList.size();s++){
		    if(proList[p2].InterfaceList[s].name==i2name)
			i2=proList[p2].InterfaceList[s].globalIndex;
		}
		proList[p1].InterfaceList[j].interfacePartnerIndex.push_back(i2);//add this to the interfaceList of partners
	    }
	}
    }//found index of each interface's partner interface.

    /*Fill in the proteinClass, now use that to fill in wholep, p_home, and Speclist.
     */
  cout <<" ------------------- "<<endl;
  cout <<" TOTAL PROTEINS: "<<proList.size()<<" TOTAL INTERFACES: "<<nTotalInterfaces<<" MAX PARNTER PER INTERFACE: "<<maxInterfacePartners<<endl;
  cout <<" ---------------------" <<endl;
  cout <<" Total IIN Edges for these proteins: "<<nTotalEdges<<". Skipping a total of: "<<NedgeRead-nTotalEdges<<" edges from the file "<<argv[3]<<" that did not have proteins listed in the file "<<argv[2]<<endl;
  cout <<"----------------------------------"<<endl;
  cout <<" PRINT OUT PROTEIN INTERFACES and PROPERTIES : "<<endl;
  cout <<"----------------------------------"<<endl;
  /*Print out interaction network and Protein properties*/
  for(p1=0;p1<proList.size();p1++){
      cout <<" Protein index: "<<p1<<" Name: "<<proList[p1].name<<" copies: "<<proList[p1].copies<<" Total interfaces: "<<proList[p1].InterfaceList.size()<<" Name, index: "<<'\t';
      for(j=0;j<proList[p1].InterfaceList.size();j++){
	  i1=proList[p1].InterfaceList[j].globalIndex;
	  cout <<proList[p1].InterfaceList[j].name<<","<<proList[p1].InterfaceList[j].globalIndex<<'\t';
      }
      cout <<endl;
      /*Below commented out, printing all interface binding partners.*/
      /*
      for(j=0;j<proList[p1].InterfaceList.size();j++){
	  cout <<proList[p1].name<<'\t'<<proList[p1].InterfaceList[j].name<<" Partners: "<<proList[p1].InterfaceList[j].interfacePartnerName.size()<<endl;
	  for(int k=0;k<proList[p1].InterfaceList[j].interfacePartnerName.size();k++){
	       cout <<proList[p1].InterfaceList[j].interfacePartnerName[k]<<","<<proList[p1].InterfaceList[j].interfacePartnerIndex[k]<<'\t';
	  }//interface partners
	  cout <<endl;
      }//interfaces
      */
  }//proteins

  string *iNames=new string[nTotalInterfaces];
  int **Speclist=new int*[nTotalInterfaces];
  double **Stoichlist=new double*[nTotalInterfaces];
  for(i=0;i<nTotalInterfaces;i++){
      Speclist[i]=new int[maxInterfacePartners];
      Stoichlist[i]=new double[maxInterfacePartners];
  }
    int *numpartners=new int[nTotalInterfaces];
    Protein *wholep=new Protein[proList.size()];
    /*assign a unique interface to each protein*/
    int Nif=nTotalInterfaces;
    
    int *p_home=new int[Nif];//this reverses and tells you what protein a given interface belongs to
    t=0;
    /*Loop over all proteins, figure out the indexes of each interface interaction using the global indices.*/

    for(p1=0;p1<proList.size();p1++){
	int nint=proList[p1].InterfaceList.size();
	wholep[p1].ninterface=proList[p1].InterfaceList.size();
	for(j=0;j<nint;j++){
	    i1=proList[p1].InterfaceList[j].globalIndex;
	    wholep[p1].valiface[j]=i1;//interfaces on this protein.
	    iNames[i1]=proList[p1].name;
	    iNames[i1].append(".");
	    iNames[i1].append(proList[p1].InterfaceList[j].name);
	    p_home[i1]=p1;//home protein for this interface.
	    int numPartner=proList[p1].InterfaceList[j].interfacePartnerName.size();
	    numpartners[i1]=numPartner;
	    for(int k=0;k<numPartner;k++){//umpartners[i1];k++){
		int i2=proList[p1].InterfaceList[j].interfacePartnerIndex[k];
		Speclist[i1][k]=i2;
		Stoichlist[i1][k]=proList[p1].InterfaceList[j].stoichiometry[k];
	    }
	}

	

    }
    cout <<"----------------------------------"<<endl;
    cout <<" PRINT OUT INTERFACE INTERACTION NETWORK "<<endl;
    cout <<"----------------------------------"<<endl;
    /*Now write out Speclist*/
    for(i=0;i<nTotalInterfaces;i++){
	cout <<i<<" Name: "<<iNames[i]<<" Npartner: "<<numpartners[i]<<'\t'<<" Partners: ";
	for(j=0;j<numpartners[i];j++){
	    //cout <<Speclist[i][j]<<'\t';
	    int id=Speclist[i][j];
	    cout<<iNames[id]<<'\t';
	}
	cout <<endl;
	
    }
    cout <<"----------------------------------"<<endl;
    
    vector<double> bin;
    //double *cmf = new double[MAXBINS];
    vector<double> cdf;
    
  t=1;
  bin.push_back(0);//first value of copy number binds
  cdf.push_back(0);//first value of CDF
  cout <<"---------------"<<endl;
  int nruns=plist.nruns_rand; //Number of random runs
  if(argc==5){  
      
      ifstream distfile(argv[4]);
      double btmp;
      double ctmp;
      /*Read in conc distribution*/
      cout <<"READING IN CELL TYPE COPY NUMBER DISTRIBUTION FROM THE FILE: "<<argv[4]<<endl;  
      distfile.ignore(400,'\n');//ignore column name titles.
      while(distfile >> btmp >> ctmp){
	  bin.push_back(btmp);
	  cdf.push_back(ctmp);
	  t++;
      }
      distfile.close();
      cout <<"TOTAL LINES in Copy Number Distribution for Cell Type: "<<t<<endl;
      // for(t=0;t<cdf.size();t++){
// 	  cout <<bin[t]<<' '<<cdf[t]<<endl;
//       }

      /*make sure max_complex is big enough to accomodate the large protein copy numbers in this distribution.*/
      if(plist.max_complex < bin[bin.size()-1] && plist.max_complex>0){
	  cout<<"---------------------WARNING-----------------------"<<endl;
	  cout <<" WARNING: max_complex value is too low--limits the complexes and they might not be able to sum up to maximum copy numbers observed in this cell type. "<<'\t';
	  cout <<"setting max_complex to the largest copy number in the cell-type: ";
	  plist.max_complex=bin[bin.size()-1];
	  cout<<plist.max_complex<<endl;
	  cout<<"--------------------------------------------"<<endl;
      }
  }else{
      nruns=0;
      cout<<" NO CELL TYPE DISTRIBUTION FILE READ IN. NOT PERFORMING STOICHIOMETRIC BALANCE FOR ANY RANDOMIZED COPY NUMBERS--NO P-VALUE CALCULATION."<<endl;
  }
     
  /*Create the matrix A*/
  int Ncomplex;
  string *edgenames=new string[nTotalEdges];
  Ncomplex=0; //Need to account for extra edges in IIN
  for(i=0;i<Nif;i++){
    for(j=0;j<numpartners[i];j++){
      if(Speclist[i][j]>=i){//Avoid counting an edge twice
	
	  edgenames[Ncomplex]=iNames[i];
	  edgenames[Ncomplex].append("::");
	  edgenames[Ncomplex].append(iNames[Speclist[i][j]]);
	  Ncomplex++;    
      }
    }
  }
  cout <<"-------------------"<<endl;
  cout << "IIN edges: " << Ncomplex << endl;
  if(Ncomplex != nTotalEdges){
      cerr <<" ERROR, EXITING: TOTAL DEFINED EDGES IN SPECLIST DOES NOT MATCH NUMBER DEFINED ABOVE "<<endl;
      cerr<<" Nc: "<<Ncomplex<<" doesn't match: "<<nTotalEdges<<endl;
      exit(1);
  }
  double *A=new double[Nif*Ncomplex];
  // cout <<"matrix size: "<<Nif*Ncomplex<<endl;
  for(i=0;i<Nif*Ncomplex;i++)
    A[i]=0;
    
  build_Amatrix_with_Stoich(Nif, A, numpartners, Speclist, Stoichlist);
  
  /*A[i*Nif+j] i is column, j is the row, each row is a different interface, each column a different complex.*/
  /* cout <<"-------------------"<<endl;
  cout <<" Write A MATRIX TRANSPOSE (First row is all interfaces in complex 1) to amatrixT.out file: "<<endl;
  ofstream amatrix("amatrixT.out");
  for(i=0;i<Nif*Ncomplex;i++){
      if((i+1)%Nif==0)
	  amatrix <<A[i]<<endl;
      else
	  amatrix <<A[i]<<'\t';
  }
  cout <<"-------------------"<<endl; 
  */
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
    cupp[i]=plist.max_complex;//In reality , this could be set for each element of the array, instead of using same value for all!
    if(plist.max_complex>0)
	icupp[i]=1;
    else
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
  
  for(i=0;i<Nconstrain;i++){
      for(j=0;j<wholep[constrain[i]].ninterface;j++){
	  ntmp=wholep[constrain[i]].valiface[j];
	  Z[ntmp*Nif+ntmp]=1;
      }
  }
  
  /*
  for(i=0;i<Npro;i++)
      constrain[i]=i;
  //Set Z to identity matrix
  for(i=0;i<Nif;i++){      
      Z[i*Nif+i]=1;
  }
  */
  
  

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
  if(nruns>0){
      cout <<"--------------"<<endl;
      if(plist.flagread==1)
	  cout<<"For p-value: Using Copy number sampled from the distribution in: "<<argv[4]<<endl;
      else
	  cout<<"For p-value: Using Copy Numberss shuffled from the real CNs in file: "<<argv[2]<<endl;
  }
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
  solvefile << "Protein\tInterface\tR\tR'\tR''\tvar for R'"<<endl;

  /***** Part 3 *****/
  
  /*Now that everything has been defined, do qp_solve*/
  /*We want to do this both with the real and random protein concentrations
    The first run (nc=0) will do the real, observed copy numbers for the proteins.
   */

  double pJSD=0;
  double pCSD=0;
  double realCSD,realJSD;
  
  int check=0;//Use this if some proteins not constrained
  double *rval = new double[Npro];
  //Set rvals to abund for first run
  for(i=0;i<Npro;i++)
    rval[i]=abund[i];
  cout<<"----------"<<endl;
  cout <<"NRUNS RANDOM COPY NUMBERS, to evaluate p-value: "<<nruns<<endl;

  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------------"<<endl;
  cout<<"---------------------PERFORM STOICHIOMETRIC BALANCE-------------------"<<endl;
  for(int nc=0;nc<nruns+1;nc++){
      start_timer(&qptime);
      
    if(nc>0){//Use random copy numbers
	cout <<"-------------------"<<endl;
	cout <<"Solving balance for the RANDOMIZED COPY NUMBERS,  to define a p-value. Run Number: "<<nc<<endl;
	if(plist.flagread==1)//Substituting flagread to determine what randomization method to use
	    sample_dist(Npro, bin,  cdf, rval);
	else
	    //Should only use shuffle when not constraining proteins!
	    shuffle_abund(Npro,abund,rval,Nconstrain,constrain);
	/*cout <<"INITIAL COPY NUMBERS : "<<endl;
	//Set "abund" to rvals
	for(i=0;i<Npro;i++){
	    cout <<rval[i]<<endl;
	}//abund[i]=1.0*rval[i];
	*/
      //}
    }else{
	cout <<"-------------------"<<endl;
	cout <<"--------SOLVING BALANCE FOR THE REAL, OBSERVED PROTEIN COPY NUMBERS------"<<endl;
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

    /*First run is the Real Observed Copy numbers, write out full solution*/

    
    if(nc==0){
	/*for(i=0;i<Nif;i++){
	solvefile << p_home[i] << "\t" << R[i] << "\t"<< Rp[i] << "\t" << Rpp[i] << "\t" << ivars[p_home[i]] << endl;
      }
	*/
	for(p1=0;p1<proList.size();p1++){
	    for(i1=0;i1<proList[p1].InterfaceList.size();i1++){
		int index=proList[p1].InterfaceList[i1].globalIndex;
		solvefile<<proList[p1].name<<"\t"<<proList[p1].InterfaceList[i1].name<<"\t"<<R[index]<<'\t'<<Rp[index]<<'\t'<<Rpp[index]<<'\t'<<ivars[p1]<<endl;
	    }//all interfaces on pro p1
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

void sample_dist(int npro, vector<double> bin, vector<double> cdf, double *rval)
{
  /*This function samples protein distributions from an input file 
    CellType_CopyNumberDistribution.txt*/


//   const gsl_rng_type * T;
//   gsl_rng * r;
//   gsl_rng_env_setup();
//   T = gsl_rng_default;
//   r = gsl_rng_alloc(T);

//   timeval tim;
//   gettimeofday(&tim,0);
//   double t1=tim.tv_sec + tim.tv_usec;
//   int seed=t1;
//   gsl_rng_set(r, t1);
  
  //int npro = atoi(argv[1]);
  
  int i;
  //int *rval = new int[npro];
  double rnum;
  int t;
  for(i=0;i<npro;i++){

    t=0;
    rnum=1.0*rand()/(1.0*RAND_MAX);//gsl_rng_uniform(r);
    //cout << "rnum is: " << rnum << '\t';//endl;
    while(rnum>cdf[t])
      t++;

    t--;
    
    double fract=(cdf[t+1]-rnum)/(cdf[t+1]-cdf[t]);
    rval[i]=int(fract*(bin[t+1]-bin[t])+bin[t]);
    //cout <<" fract: "<<fract<<" t (index): "<<t<<" rval: "<<rval[i]<<" Bin: "<<bin[t]<<endl;
  }
  //infile.close();
  
  //ofstream outfile;
  //outfile.open("protnums_rand.inp",ios::out|ios::trunc);
  //for(i=0;i<npro;i++)
  //outfile << rval[i] << endl;

  //outfile.close();

  //  gsl_rng_free(r);
  
}
void shuffle_abund(int npro, double *abund, double *rval, int Nconstrain, int *constrain)
{
  //Shuffle the copy number values for each protein.

  //Set random generator seed
    /*
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
    */
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
      rnum=1.0*Nconstrain*rand()/(1.0*RAND_MAX);//gsl_rng_uniform(r);
    j=int(rnum);
    //Swap a[i2] and a[j2];
    j2=constrain[j];
    i2=constrain[i];
    temp = rval[i2];
    rval[i2]=rval[j2];
    rval[j2]=temp;

  }

  //  gsl_rng_free(r);

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
      //cout << i1 << ": " << R[i1] << endl;
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

