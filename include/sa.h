#ifndef sa_h
#define sa_h
#include "atom.h"
#include <iostream>
namespace saconst{
extern double sa_temp;
extern double sa_ratio;
extern double sa_max;
extern double sa_eweight;
extern double sa_fweight;
extern double sa_sweight;
extern double sa_nt;
extern double sa_ns;
extern int sa_refindex;
extern int sa_atom_num;
}
void SimulatedAnnealing(double (*PenaltyFunc)(double*, box*,int,int,int),
			box* system,
            double *xacc,
			int    N,
			int    NT,
			int    NS,
			int    sa_max,
			double sa_temp,
			double sa_ratio,
			double *vm,
			double *ub,
			double *lb,
			double *c
);

#endif
