

void calc_ratio(int nwhole, ppidata *ppi, double *ratio);
void order_ratio(int nwhole, double *ratio, long unsigned int *index);
void calc_ratio_mc(int ninterface, int **Speclist, int *numpartners, int *nbor);
double calc_avg_ratio(int ninterface, int **Speclist, int *numpartners);
