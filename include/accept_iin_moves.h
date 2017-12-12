
void accept_swap(double newfit, double &oldfit, int ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist);
void accept_swap_iface(double newfit, double &oldfit, int ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist, int nwhole,  Protein *wholep, Protein *wholetemp, int *p_home );
void accept_ifacechg(double newfit, double &oldfit, int ntmpinterface, int &ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist, int nwhole, Protein *wholep, Protein *wholetemp);
//void accept_ifacechg(double newfit, int ntmpinterface, int &ninterface, int *numpartners, int **Speclist, int *tmppartners, int **templist, int nwhole, Protein *wholep, Protein *wholetemp);
