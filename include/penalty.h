#ifndef penalty_h
#define penalty_h
#include "atom.h"
#include "readpara.h"

void indexbvvmap(int** bvvmatrixmap,int tick,int& i,int& j);
void indexchargemap(int* chargemap,int tick,int& i);
void map(double* xp);
void mapToXp(double* xp);
double PenaltyFunc(double* xp, box* system,int number,int index,int databasetick);
#endif
