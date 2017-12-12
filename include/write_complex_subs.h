
void write_status2(ofstream &outfile, gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, int iter, int Nproonly);
void complex_sizes(ofstream &outfile, ofstream &outfile1, gillComplex *ind_com, gillFullmol *bases, Parms &plist, int *p_home, int iter, double currtime, int Nproonly, int *Ncopy);
void complex_sizes_skippip(ofstream &outfile, ofstream &outfile1, gillComplex *ind_com, gillFullmol *bases, Parms &plist, int *p_home, int iter, double currtime, int Nproonly, int *Ncopy);
void write_propartnernums(ofstream &outfile, gillProtein *wholep, gillFullmol *bases, Parms &plist, int *Ncopy, int *Nsum,  int **Rlist, int **myrxn, int *p_home, int iter, double currtime);
