#ifndef image_h
#define image_h
#include <math.h>
/*input is the atom array and size is the size of how many atoms*/
/*I use this code to construct the virtual image of this simulation unit cell,input is the simulating atoms in the box. size if how many atoms are in the cell. p is the periodical boundary condition, cutoff is the cutoff of this type of force*/
atom* imageall(atom* input,int size,double* p,double cutoff,int& virt);
/*shift the simulation box using shiftv egg{1,0,0} means shift right one period*/
void shift(atom* input,int size,double* p,int* shiftv);
#endif
