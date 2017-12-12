#include "pro_classes.h"
#include "read_yeast_conc.h"


void read_yeast_conc(double *cumul, double *abund, int &length, double &conc_mean)
{
  int N2=6234;
  char fname[200];

  sprintf(fname, "yeast_Ghaemm_exp.txt");

  ifstream infile2;
  infile2.open(fname);
  if(!infile2){
    cerr<<"Can't find protein abundance datafile: "<<fname<<endl;
    exit(1);
  }    
    
  infile2.ignore(6000,'\n');
  double bsize=100;
  //  int nbins=16000;
  string skip;
  string *cname=new string[N2];
  double *conc=new double[N2];
  
  double *hist=new double[NBINS];
  //double *cumul=new double[nbins];
  //double *abund=new double[nbins];
  double Ntotal=0.0;
  int id;
  char *str=new char[100];
  char *str1=new char[100];
  char *str2=new char[100];
  char *str3=new char[100];
  int i, k;
  conc_mean=0.0;
  for(i=0;i<NBINS;i++)
    hist[i]=0;
  ofstream nfile("yeast_conc.out");
  for(k=0;k<N2;k++){
    infile2 >>cname[k]>>skip>>str>>str1>> str2>>str3;
    conc[k]=-1;
    if(str[0]!='-')
      conc[k]=atof(str);
    //cout <<cname[k]<<' '<<skip<<' '<<conc[k]<<' '<<str1<<' '<<str2<<' '<<str3<<endl;
    /*generate a histogram of the concentrations*/
    if(conc[k]>0){
      id=int(conc[k]/(1.0*bsize));
      hist[id]++;
      Ntotal+=1.0;
      nfile<<cname[k]<<' '<<conc[k]<<endl;
      conc_mean+=conc[k];
    }
  }
  conc_mean/=(1.0*Ntotal);
  cout <<"read in concentration data. Mean concentration:  "<<conc_mean<<endl;
  cout <<"calculate distribution "<<Ntotal<<endl;
  ofstream pfile("Probs_abundance.out");
  pfile.precision(12);
  double sum=0.0;
  int s=1;
  cumul[0]=0;
  abund[0]=0;
  for(i=0;i<NBINS;i++){
    hist[i]/=(Ntotal*1.0*bsize);
    sum+=hist[i]*bsize;

    if(hist[i]>0){
      cumul[s]=sum;
      abund[s]=(i+0.5)*bsize;
      pfile<<(i+0.5)*bsize<<'\t'<<hist[i]<<'\t'<<cumul[s]<<endl;
      s++;
    }
  }
  int Nprobs=s;
  cout <<"length of probability vector: "<<Nprobs<<endl;
  length=Nprobs;
  
}
