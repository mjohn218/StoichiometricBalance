
void write_speclist(int ninterface, int *numpartners, int **Speclist );
void writetofile_speclist(ofstream &filename, int ninterface, int *numpartners, int **Speclist );
void write_ppi(int nwhole, Protein *wholep);
void writetofile_ppi(ofstream &filename, int nwhole, Protein *wholep);
void write_Edgemat(int nwhole, int **Edgemat, ppidata *ppi, ofstream &efile);

