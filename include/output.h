#include "atom.h"
#include <string>
#include <iostream>
#include "readpara.h"
#include <fstream>
#include "sa.h"
#include <cmath>
#include <iomanip>
void writeoutput(box* input,int size,int tick,int saiter,std::string deoptfile,std::string dfoptfile,std::string dsoptfile);
void write_opt_parameter(std::fstream& fs);
void write_defile(int databasetick);
