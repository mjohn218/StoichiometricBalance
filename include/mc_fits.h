double mc_fit(int ninterface, int *nbor, int *numpartners, double beta);
double mc_fit2(int ninterface, int *nbor, int *numpartners, double beta, int nwhole, Protein *wholep, double mu);
double mc_fit3(int ninterface, int *nbor, int *numpartners, double beta, int nwhole, Protein *wholep, double mu, double gridcoeff, double ksi, double edgediff);
double mc_fit4(int ninterface, double *grid_cofs, double *clocals, double ksi, double beta, int nwhole, Protein *wholep, double mu, int *numpartners, double edgediff);
double mc_fit5(int ninterface, double *grid_cofs, double *clocals, double ksi, double beta, int nwhole, Protein *wholep, double mu, int *numpartners, double edgediff);
double mc_fit6(int ninterface, double *grid_cofs, double *clocals, double ksi, double beta, int nwhole, Protein *wholep, double mu, int *numpartners, double edgediff, int **Speclist, double zeta);
