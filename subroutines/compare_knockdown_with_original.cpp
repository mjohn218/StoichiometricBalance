#include <vector>
#include "ProteinClass.h"
#include "pro_classes.h"

#include "constrainParms.h"
#include "construct_net.h"
#include "calc_abund_target_metrics.h"

void compare_knockdown_with_original(ofstream &deltaFile, ofstream &matrixFile, ofstream &matrixFileRB, string label, vector<ProteinClass> proListKD, vector<int> proKDID, vector<ProteinClass> proList, vector<double> realKD, vector<double> iFaceBalancedKD, vector<double> proBalancedKD, vector<double> iFaceBalanced, vector<double> proBalanced)
{
  
  /*Because the Original has more proteins in its list, need to remove those entries from its proBalanced List.*/
  int p1, i1;
  int Nif=0;
  vector<double> origReal;
  vector<double> origIBalanced;
  vector<double> origProBalanced;
  vector<double> localReal;
  vector<double> localIB;
  vector<double> localPB;

  vector<double> origRBratio;
  vector<double> KDRBratio;
  vector<double> balancePerPro;
  double rb;
  /*Need to find the matching protein/interfaces in both lists, they are both out-of-order*/
  int currI=0;
  for(p1=0;p1<proList.size();p1++){
    /*if p1 is in the KD list, need to skip it.*/
    bool keep=true;
    for(int s=0;s<proKDID.size();s++){
      if(proKDID[s]==p1){
	cout <<"exclude protein: "<<p1<<endl;
	keep = false;
	break;
      }
    }
    /*if p1 is not in the list for proListKD, due to lost edges, for example, then skip it.*/
    bool found = false;
    for(int p2=0;p2<proListKD.size();p2++){
      if(proListKD[p2].name==proList[p1].name){
	found = true;
	cout <<" Both lists contain: "<<proListKD[p2].name<<endl;
	break;
      }
    }
    int index;
    int index2;
    bool matched;
    double c;
    if(keep==true && found ==true){
      //add these to the list of copy numbers.
      for(i1=0;i1<proList[p1].InterfaceList.size();i1++){
	index=proList[p1].InterfaceList[i1].globalIndex;
	/*Find the matching index in the KD list*/
        matched=false;
	for(int p2=0;p2<proListKD.size();p2++){
	  if(proListKD[p2].name==proList[p1].name){
	    //we found the matching protein, now find the matching interface.
	    //cout <<"found matching proteins at indexes: "<<p1<<" "<<p2<< ' '<<proListKD[p2].name<<" copies: "<<proList[p1].copies<<' '<<proListKD[p2].copies<<endl;
	    for(int i2=0;i2<proListKD[p2].InterfaceList.size();i2++){
	      if(proListKD[p2].InterfaceList[i2].name==proList[p1].InterfaceList[i1].name){
		matched=true;
		origReal.push_back(proList[p1].copies);//this is to verify that they are in the right order!
		origIBalanced.push_back(iFaceBalanced[index]);
		origProBalanced.push_back(proBalanced[index]);
		
		index2=proListKD[p2].InterfaceList[i2].globalIndex;
		//cout <<"Found matching interfaces, "<<index2<<" proBalanced copies:" <<proBalancedKD[index2]<<" origCopies: "<<proBalanced[index]<<" name: "<<proListKD[p2].InterfaceList[i2].name<<endl;
		double a=proListKD[p2].copies;
		double b=iFaceBalancedKD[index2];
		c=proBalancedKD[index2];
		localReal.push_back(a);
		localIB.push_back(b);
		localPB.push_back(c);
		//cout <<" I: "<<Nif<<" protein: "<<proList[p1].name<<" CN "<<origReal[Nif]<<" Ibalance "<<origIBalanced[Nif]<<" probalance: "<<origProBalanced[Nif]<<endl;
		//cout <<" withKD "<<proListKD[p2].name<<" KDCopies "<<localReal[Nif]<<" KDIbalance: "<<localIB[Nif]<<" KDprobalance: "<<localPB[Nif]<<endl;
		break;
	      }
	    }
	    break;
	  }
	}//done finding index of KD protein
	if(matched==true)Nif++;
      }
      /*Store the proBalanced in an array of size NproKD*/
      balancePerPro.push_back(c);
      
      //take ratio of the R/B.
      //exclude any proteins with copies<0 (no known values)
      if(proList[p1].copies>0){
	rb=1.0*proList[p1].copies/(1.0*proBalanced[index]);
	origRBratio.push_back(rb);//This ratio only includes proteins with known Real copies.
	
	//cout <<"orig RBratio: "<<p1<<' '<<proList[p1].name<<' '<<rb<<" copie: "<<proList[p1].copies<<" balanced: "<<proBalanced[index]<<" KD balanced: "<<balancePerPro[currI]<<endl;
      }
      currI++;
    }//skip proteins that were KDed or remoeved from the KD lists.
  }//done re-constructing the balanced lists from original-KD.
  cout <<" MATCHED Ninterfaces: "<<Nif<<endl;
  if(balancePerPro.size() != proListKD.size()){
    cerr <<" ERROR IN compare_knockdown_with_original. MISMATCHED NUMBER OF PROTEINS R/B RATIOS CALCULATED VS NUMBER OF TOTAL PROTEINS IN PROLISTKD "<<endl;
    exit(1);
  }
  /*Construct RB ratio for the KD*/
  int nProRB=0;
  int index;
  for(p1=0;p1<proListKD.size();p1++){
    index=proListKD[p1].InterfaceList[0].globalIndex;
    if(proListKD[p1].copies>0){
      
      rb=proListKD[p1].copies/proBalancedKD[index];
      KDRBratio.push_back(rb);
      //cout <<"KD RBratio: "<<p1<<' '<<proListKD[p1].name<<' '<<rb<<endl;
      nProRB++;
    }
  }

  /*print out*/
  //  cout <<" Length of original RB ratio: "<<origRBratio.size()<<" Length of KD: "<<KDRBratio.size()<<endl;
  //for(int i=0;i<KDRBratio.size();i++){
  //cout<<origRBratio[i]<<' '<<KDRBratio[i]<<endl;
  //}
  //cout <<" Balanced copy numbers "<<endl;

  double CSD[3];
  CSD[0]=calc_chisquare(Nif,origReal,localReal);//THIS SHOULD BE ZERO, IDENTICAL COPIES!!
  if(CSD[0]!=0){
    cout <<" WARNING: ERROR!-------------- "<<endl;
    cout <<" The comparisons between the copy numbers original and after knockdown are not in the right order! "<<endl;
    deltaFile<<" WARNING: ERROR ON THIS DISTRIBUTION--OUT OF ORDER ENTRIED "<<endl;
  }
  /*When comparing real observed copies with balanced, ignore entries where real=-1. */
  vector<double> localTemp;
  for(int i=0;i<Nif;i++){
    if(localReal[i]<0)
      localTemp.push_back(-1);
    else
      localTemp.push_back(localPB[i]);
    //    cout <<" Interface: "<<i<<' '<<localReal[i]<<' '<<localTemp[i]<<endl;
  }
  
  CSD[0]=calc_chisquare(Nif, localReal, localTemp);//this should be the same as is measured in solve_balanced_from_inputs
  CSD[1]=calc_chisquare(Nif,origIBalanced,localIB);
  CSD[2]=calc_chisquare(Nif,origProBalanced,localPB);
  
  double ENTR[3];
  //Calculating Jensen Shannon distance, rather than KL divergence
  ENTR[0]=calc_dist_jensenshannon(Nif,localReal, localPB);
  ENTR[1]=calc_dist_jensenshannon(Nif,origIBalanced, localIB);
  ENTR[2]=calc_dist_jensenshannon(Nif,origProBalanced, localPB);
  
  /*Compare distances between the R/B ratios*/
  double chiRB=calc_chisquare(KDRBratio.size(), origRBratio, KDRBratio);
  double jsdRB=calc_dist_jensenshannon(KDRBratio.size(), origRBratio, KDRBratio);

  //Write out entries to file
  deltaFile << label << "\t" <<CSD[0]<<"\t"<<ENTR[0]<<"\t"<< CSD[1] << "\t" << CSD[2] << "\t" << ENTR[1] << "\t" << ENTR[2] << "\t"<<chiRB<<"\t"<<jsdRB<<endl;
  
  /*Write out full solutions to the per Protein files.*/
  int np=0;
  matrixFile<<label<<'\t';
  matrixFileRB<<label<<'\t';
  for(p1=0;p1<proList.size();p1++){
    /*if p1 is in the KD list, need to skip it.*/
    bool keep=true;
    for(int s=0;s<proKDID.size();s++){
      if(proKDID[s]==p1){
	cout <<"exclude protein: "<<p1<<endl;
	keep = false;
	break;
      }
    }//checking if we have this protein in our list
    /*if p1 is not in the list for proListKD, due to lost edges, for example, then skip it.*/
    bool found = false;
    for(int p2=0;p2<proListKD.size();p2++){
      if(proListKD[p2].name==proList[p1].name){
	found = true;
	cout <<" Both lists contain: "<<proListKD[p2].name<<endl;
	break;
      }
    }
    
    if(keep==true && found==true){
      //print out the Balanced copies, they need to be separately indexed, since the vector is length NproKD
      matrixFile<<balancePerPro[np]<<'\t';
      if(proList[p1].copies>0){
	//print out the RB ratio
	matrixFileRB<<proList[p1].copies/balancePerPro[np]<<'\t';
      }
      np++;
    }else{
      matrixFile<<-1<<'\t';
      if(proList[p1].copies>0)
	matrixFileRB<<-1<<'\t';
    }
  }
  matrixFile<<endl;
  matrixFileRB<<endl;
  if(np != proListKD.size()){
    /*Check--should have only printed out RB ratios for the length of proListKD */
    cerr <<" ERROR IN compare_knockdown_with_original. MISMATCHED NUMBER OF PROTEINS R/B RATIOS CALCULATED VS NUMBER OF TOTAL PROTEINS IN PROLISTKD "<<endl;
    exit(1);
  }
  
}
