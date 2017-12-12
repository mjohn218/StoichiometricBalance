

void name_interfaces(int Nifaces, string *interface_names, int *phome, int *ihome, string *pronames);
void name_products(int Nrxn, int **Rlist, int *rxntype, int Nifaces, ProductName *complexnames, string *proname, int *phome, int *ihome);
void write_conc_names(ofstream &corrfile1, ofstream &corrfile2, int Nifaces, int Nspecies, string *interface_names, ProductName *complexnames, ofstream &keyfile);
void write_conc_names_nonzero(ofstream &corrfile1, ofstream &corrfile2, int Nifaces, int Nspecies, string *interface_names, ProductName *complexnames, ofstream &keyfile);
