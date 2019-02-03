#ifndef interface_h
#define interface_h
#include <string>
#include "atom.h"
double** getparameter(std::string pair_style,std::string file);
atom* configuration(std::string file);
#endif
