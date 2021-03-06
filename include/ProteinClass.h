#include <string>
#include <vector>
using namespace std;

class Interface{

public:
    string name;
    //string proteinHome;
    int globalIndex;
    int nPartners;
    //    vector<string> partnerList;
    vector<string> proPartnerName;
    vector<string> interfacePartnerName;
    vector<int> proPartnerIndex;
    vector<int> interfacePartnerIndex;//global index
    Interface(int index);//constructor
    vector<double> stoichiometry;
};


class ProteinClass{
public:
    int nInterfaces;
    string name;
    //    vector<string> iFaceNames;
    vector<Interface> InterfaceList;//has properties of the above class
    int nProPartners;
    //vector<string> proPartnerNames;
    double copies;
    bool isConstrained;
    int index;
    ProteinClass(double n_copies);//constructor
    ProteinClass();//default constructor
};
