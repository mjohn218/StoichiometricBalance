
void construct_network_from_inputs(int NedgeRead, int &nTotalInterfaces, int &maxInterfacePartners, int &nTotalEdges, vector<ProteinClass> &proList, vector<string> &edge1Interface, vector<string> &edge2Interface, vector<string> &edge1Protein, vector<string> &edge2Protein, vector<string> &stoichiometry, vector<string> filenames);
void construct_formatted_interface_arrays(int nTotalInterfaces, vector<ProteinClass> proList, Protein *wholep, int *numpartners, string *iNames, int **Speclist, double **Stoichlist, int *p_home);
void define_edgenames(int &Ncomplex, int Nif, int *numpartners, int **Speclist, string *edgenames, int nTotalEdges, string *iNames);
void create_new_prolist( vector<ProteinClass> &proList, vector<ProteinClass> &proListKD, vector<int> proKDID, int nProKD);
void reset_arrays(vector<string> &knockDownNames, vector<int> &proKDID, vector<ProteinClass> &proListKD, vector<double> &realKD, vector<double>&iFaceBalancedKD, vector<double> &proBalancedKD, vector<double> &stdBalancedKD);
void compare_knockdown_with_original(ofstream &deltaFile, ofstream &matrixFile, ofstream &matrixFileRB, string label, vector<ProteinClass> proListKD, vector<int> proKDID, vector<ProteinClass> proList, vector<double> realKD, vector<double> iFaceBalancedKD, vector<double> proBalancedKD, vector<double> iFaceBalanced, vector<double> proBalanced);

void init_matrix_file(ofstream &matrixFile, ofstream &matrixFileRB, vector<ProteinClass> proList, vector<double> proBalanced);
