#include <fstream>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cstdlib>

int trand();
int mutate_connections(int nwhole, int *numpartners, int **Speclist, Protein *wholep);
int mutate_edge(int nwhole, int *numpartners, int **Speclist, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag, int *p_home);
int mutate_interfaces_rev(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag);
int combine_interfaces_rev(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, double &pf, double &pb);
int split_interfaces_rev(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, double &pf, double &pb);

int split_type(int nedge);
int split_type_self(int nedge, int &flagsep);

int add_edge(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int maxne, int *selfflag, double &pf, double &pb, int PPIedge);
int delete_edge(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, int PPIedge);

int find_copies(int nwhole, int N, int *numpartners, Protein *wholep, int **Speclist, int *Adj, int i1, int *p_home, int bait, int &copy);

int mutate_connections_orig(int nwhole, int *numpartners, int **Speclist, Protein *wholep);
int mutate_edge_orig(int nwhole, int *numpartners, int **Speclist, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag, int *p_home);
int mutate_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *selfflag);
int combine_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, double &pf, double &pb);
int split_interfaces_rev_orig(int nwhole, int ninterfaces, int *numpartners, int **Speclist, int *Adj, Protein *wholep, ppidata *ppi, double &pgen_ratio, int *p_home, int maxni, int *selfflag, double &pf, double &pb);
