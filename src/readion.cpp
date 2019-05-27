#include "readion.h"
#include "atom.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include "readpara.h"
#include "penalty.h"
#include <mpi.h>
#include <math.h>
/*number is the number of atoms,boxnumber is the how many number of structures are there, ref is the reference structure. cutoff is the cut-off the iteractive force.*/
box* readion(std::string inputfile,int databasetick,int& boxnumber,int& ref,double cutoff){
	std::fstream fs;
	fs.open(inputfile.c_str(),std::fstream::in);
	std::string line;
	std::istringstream stream1;
	getline(fs,line);
	stream1.str(line);
	int flag=0;
	stream1>>flag;
	boxnumber=flag;
	stream1.clear();
	getline(fs,line);
	stream1.str(line);
	stream1>>flag;
	int number=flag;
	stream1.clear();
	/*
	 * For N MPI processors, if i!=N-1, this processor should have ceil(m/N) data point.
	 * for i=N-1, this processor should have min(ceil(m/N),m-(N-1)*ceil(m/N)) data point.
	 * when garthering information 
	 * */
	int world_rank,mpi_size;
	control::dftenergy[databasetick]=new double [boxnumber];
	control::mdenergy[databasetick]=new double [boxnumber];
	control::diffenergy[databasetick]=new double [boxnumber];
	MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
	int box_ave=floor((boxnumber+0.0000001)/mpi_size);
	int remain=boxnumber%mpi_size;
	int box_size_local;
	ref=0;
	double e_ref=1e25;
	int local_struct_ref=0;
	if(world_rank<remain){
		box_size_local=box_ave+1;
	}
	else{
		box_size_local=box_ave;
	}
	box* ionall=new box[box_size_local];
	atom* atomconfig;
	double* period=new double[3];
	double** stress_dft;
	double dftenergy;
	double weight;
	/************************set the type tick**************************/
	/*the type tick go from 0,1,2,3,4,5,6,7,.....*/
	int* type_tick=new int [number];
	for(size_t i=0;i<number;i++){
		getline(fs,line);
		for(size_t j=0;j<species::spe.size();j++){
			if(line.find(species::spe[j])!=std::string::npos){
				type_tick[i]=j;
				break;
			}
		}
	}
	if(world_rank==0){
	std::cout<<"---------------------------------------STARTING READING ION COORDINATES----------------------------"<<std::endl;
  std::cout<<"There are how many atoms "<<number<<std::endl;
	std::cout<<"you are READING "<<inputfile<<" ionic files"<<std::endl;
	std::cout<<"the input sequence for atom is: "<<std::endl;
	for(size_t i=0;i<number;i++){
		if(i%5==0){
			std::cout<<std::endl;
		}
		std::cout<<type_tick[i]+1<<"\t";
	}
	std::cout<<std::endl;
	std::cout<<"---------------------------------------------------END--------------------------------------------"<<std::endl;
	}
	/*******************************************************************/
	for(size_t tick=0;tick<boxnumber;tick++){
		getline(fs,line);
		atomconfig=new atom [number];
		for(size_t j=0;j<number;j++){
			getline(fs,line);
			stream1.clear();
			stream1.str(line);
		/*reading atomic positions*/
			for(size_t k=0;k<3;k++){
				stream1>>atomconfig[j].position[k];
			}
			stream1.clear();
			/*end reading that*/
		}
		getline(fs,line);
		/*reading periodical boudary condition*/
		for(size_t k=0;k<3;k++){
				getline(fs,line);
				stream1.str(line);
				for(size_t m=0;m<=k;m++){
						stream1>>period[k];
				}
				stream1.clear();
			}
		/*end reading that*/
		getline(fs,line);
		getline(fs,line);
		getline(fs,line);
		/*readling forces*/
		for(size_t j=0;j<number;j++){
			getline(fs,line);
			stream1.str(line);
			for(size_t k=0;k<3;k++){
				stream1>>atomconfig[j].dftforce[k];
			}
			stream1.clear();
		}
			/*end reading that*/
		getline(fs,line);
		/*reading energy*/
		getline(fs,line);
		stream1.str(line);
		stream1>>dftenergy;
		if(dftenergy<e_ref){
			ref=tick;
			e_ref=dftenergy;
		}
		stream1.clear();
		/*end reading energy*/
		getline(fs,line);
		/*reading stress tensor*/
		stress_dft=new double* [3];
		for(size_t i=0;i<3;i++){
			stress_dft[i]=new double [3];
			getline(fs,line);
			stream1.str(line);
			for(size_t j=0;j<3;j++){
				stream1>>stress_dft[i][j];
			}
			stream1.clear();
		}
		/*end reading stress tensor*/
		getline(fs,line);
		/*reading weight*/
		getline(fs,line);
		stream1.str(line);
		stream1>>weight;
		for(size_t i=0;i<number;i++){
			for(size_t j=0;j<3;j++){
				atomconfig[i].position[j]=atomconfig[i].position[j]*period[j];
			}
			atomconfig[i].type=type_tick[i];
		}
		if(tick%mpi_size==world_rank){
			ionall[local_struct_ref].init(atomconfig,number,species::spe.size(),cutoff,period,control::bvvmatrix,dftenergy,stress_dft,weight);
			local_struct_ref++;
		}
		delete [] atomconfig;
	}
	delete [] period;
	delete [] type_tick;
  return ionall;
}
