#ifndef readpara_h
#define readpara_h
/*read control.PT*/
#include <string>
#include <vector>
#include "atom.h"
#include <fstream>
#include <iostream>
namespace control{
extern double** bvvmatrix;
extern int** bvvmatrixmap;
extern double* lb;
extern double* ub;
extern double* vm;
extern double* c;
extern double* charge;
extern int* chargemap;
extern int* type;
extern int pair_num;
extern int paracount_bvv;
extern int paracount_charge;
extern std::vector<int> para_site_charge_change;
extern std::vector<double> para_site_charge;
extern std::vector<std::string> site_name;
extern std::vector<double> chemical_formula;
extern double* xop;
extern std::vector<std::vector<int> > mapXpTickToBvvTick;
extern std::vector<std::vector<int> > mapXpTickToChargeTick;
extern int lastchargetick;
extern int neutral;
extern std::vector<std::string> ionfile;
extern std::vector<box*> database;
extern std::vector<int> minienergytick;/*store the minimum energy of this Ion files*/
extern std::vector<int> ionsize;/*store the structure numbers of different files*/
extern std::vector<std::string> deopt;
extern std::vector<std::string> dfopt;
extern std::vector<std::string> dsopt;
extern std::vector<double*> dftenergy;
extern std::vector<double*> mdenergy;
extern std::vector<double*> diffenergy;
}
namespace ewaldsum{
extern double cutoff;
extern double k_cutoff;
extern double alpha;
}
namespace species{
	extern std::vector<std::string> spe;
	extern std::vector<int> nametag;
    extern std::vector<int> site;/*0 for asite, 1 for bsite, 2 for osite*/ 
	extern int** num;
}
void readPT(std::string PTfile);
void readvmmap(std::fstream &fs);
void readbound(std::string boundfile);
#endif
