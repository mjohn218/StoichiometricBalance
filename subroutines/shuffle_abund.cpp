#include "pro_classes.h"
#include "functionCalls.h"

void shuffle_abund(int npro, vector<double> abund, vector<double> &rval, int nConstrain, vector<int> constrain)
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
    //This vector has already been assigned, so we can overwrite values here this way.
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
  for(i = nConstrain-1; i>0; i--){
      rnum=1.0*nConstrain*rand()/(1.0*RAND_MAX);//gsl_rng_uniform(r);
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

