
#ifndef INIT_NET_H
#define INIT_NET_H


int init_network(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int Nedge, int **Speclist, int *numpartners );
int init_network_split(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int Nedge, int **Speclist, int *numpartners, int *e1num, int *e2num);
int init_network_both(int nwhole, ppidata *ppi, Protein *wholep, int *p_home, constrainParms &plist, int Nedge, int **Speclist, int *numpartners, int *e1num, int *e2num, double *abund);

#endif /* INIT_NET_H */
