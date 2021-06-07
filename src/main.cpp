#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include "readion.h"
#include <ctime>
#include "readpara.h"
#include "penalty.h"
#include "sa.h"
#include <mpi.h>
#include <time.h>
int main(int argc,char* argv[]){
	 MPI_Init(NULL,NULL);
	 readPT(argv[1]);
	 int size_box;
	 MPI_Barrier(MPI_COMM_WORLD);
   SimulatedAnnealing(&PenaltyFunc,
			 control::database[0],
			 control::xop,
			 control::paracount_bvv+control::independentcharge,
			 saconst::sa_nt,
			 saconst::sa_ns,
			 saconst::sa_max,
			 saconst::sa_temp,
			 saconst::sa_ratio,
			 control::vm,
			 control::ub,
			 control::lb,
			 control::c);
	 /*
	 int i=0;
	 	int world_rank;
		MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	 if(world_rank==0){
		 std::cout<<"the reference tick is: "<<control::minienergytick[i]<<std::endl;
	 }
	clock_t start=clock();
	double penaltyp;
	for(size_t k=0;k<1;k++){
		penaltyp = PenaltyFunc(control::xop,control::database[i],control::ionsize[i],control::minienergytick[i]);//Zhenbang
		MPI_Barrier(MPI_COMM_WORLD);
		if(world_rank==0){
			std::cout<<"the cycle time is: "<<k<<std::endl;
		 std::cout<<"the penalty is: "<<penaltyp<<" recycle time is"<<k<<std::endl;
		}

	}
	clock_t end=clock();
		if(world_rank==0){
		 std::cout<<"the penalty is: "<<penaltyp<<std::endl;
		 std::cout<<"the time used is: "<<(double)(end-start)/CLOCKS_PER_SEC;
		}
		*/
		MPI_Finalize();
}
