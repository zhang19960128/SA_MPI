#ifndef readion_h
#define readion_h
#include "atom.h"
#include <string>
/*number is the number of simulation atoms, boxnumber is how many boxs are created in the simulation*/
box* readion(std::string inputfile,int databasetick,int& boxnumber,int& ref,double cutoff);
#endif
