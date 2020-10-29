#include "pro_classes.h"
#include "functionCalls.h"

void sample_dist(int npro, vector<double> bin, vector<double> cdf, vector<double> &rval)
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
  double conc;
  for(i=0;i<npro;i++){

    t=0;
    rnum=1.0*rand()/(1.0*RAND_MAX);//gsl_rng_uniform(r);
    //cout << "rnum is: " << rnum << '\t';//endl;
    while(rnum>cdf[t])
      t++;

    t--;
    
    double fract=(cdf[t+1]-rnum)/(cdf[t+1]-cdf[t]);
    /*This vector has already been defined, so overwriting values here.*/
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
