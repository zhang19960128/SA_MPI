#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <iomanip>
void box::computestress(){
	for(size_t k=0;k<3;k++)
		for(size_t j=0;j<3;j++){
			stress[k][j]=0.0;
		}
	for(size_t i=0;i<size;i++){
		for(size_t k=0;k<3;k++)
			for(size_t j=0;j<3;j++){
					stress[k][j]+=allatom[i].position[k]*allatom[i].force[j];
			}
	}
	for(size_t k=0;k<3;k++){
		for(size_t j=0;j<3;j++){
//			std::cout<<stress[k][j]<<" "<<"\t";
		}
//		std::cout<<std::endl;
	}
}
