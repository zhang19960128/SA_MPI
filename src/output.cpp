#include "atom.h"
#include <string>
#include <iostream>
#include "readpara.h"
#include <fstream>
#include "sa.h"
#include <cmath>
#include <iomanip>
/*
 * size is the size_box ot input box, tick is the minimum DFT energy, deoptfile is the difference energy file, dfoptfile is the difference Force file, dsoptfile is the difference stress file
 * 
 */
void writeoutput(box* input,int size_box,int tick,int saiter,std::string deoptfile,std::string dfoptfile,std::string dsoptfile){
	std::fstream dfopt,deopt,dsopt;
	dfopt.open(dfoptfile.c_str(),std::fstream::out);
	deopt.open(deoptfile.c_str(),std::fstream::out);
	dsopt.open(dsoptfile.c_str(),std::fstream::out);
	deopt<<"#"<<saconst::sa_temp<<std::endl;
	deopt<<saiter+1<<std::endl;
	deopt<<std::setw(18)<<"ABS(E(DFT)-E(MD))\t"<<std::setw(18)<<"E(DFT)\t"<<std::setw(18)<<"E(MD)"<<std::endl;
	dfopt<<std::setw(18)<<"ABS(F(DFT)-F(MD))\t"<<std::setw(18)<<"F(DFT)\t"<<std::setw(18)<<"F(MD)"<<std::endl;
	double dftetemp;
	double mdetemp;
	for(size_t i=0;i<size_box;i++){
		mdetemp=input[i].mdenergy-input[tick].mdenergy;
		dftetemp=input[i].dftenergy-input[tick].dftenergy;
		deopt<<std::setw(18)<<"\t"<<std::fabs(dftetemp-mdetemp)<<std::setw(18)<<"\t"<<dftetemp<<std::setw(18)<<"\t"<<mdetemp<<std::endl;
		for(size_t j=0;j<(input+i)->size;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<fabs(input[i].allatom[j].force[k]-input[i].allatom[j].dftforce[k]);
			}
			}
		for(size_t j=0;j<(input+i)->size;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<input[i].allatom[j].dftforce[k];
			}
			}
		for(size_t j=0;j<(input+i)->size;j++){
			for(size_t k=0;k<3;k++){
			dfopt<<std::setw(10)<<"\t"<<input[i].allatom[j].force[k];
			}
			}
		std::cout<<std::endl;
	}
}
void write_opt_parameter(std::fstream& fs){
	fs<<"the new-optimized parameters are: "<<std::endl;
	for(size_t i=0;i<control::pair_num;i++){
	  for(size_t j=0;j<12;j++){
			fs<<std::setw(8)<<std::fixed<<control::bvvmatrix[i][j]<<"\t";
		}
		fs<<std::endl;
	}
	for(std::vector<std::string>::iterator a=species::spe.begin();a!=species::spe.end();a++){
		fs<<*a<<"\t";
	}
	fs<<std::endl;
	size_t size=species::spe.size();
	for(size_t i=0;i<size;i++){
		fs<<std::setw(8)<<std::fixed<<control::charge[i]<<"\t";
	}
	fs<<std::endl;
	fs<<"the Change-Range of Parameters are: "<<std::endl;
	int tick=0;
	for(size_t i=0;i<control::pair_num;i++){
	  for(size_t j=0;j<12;j++){
			fs<<std::setw(8)<<std::fixed<<control::bvvmatrixmap[i][j]*control::vm[tick]<<"\t";
			if(control::bvvmatrixmap[i][j]==1){
				tick++;
			}
		}
		fs<<std::endl;
	}
	for(size_t i=0;i<control::site_name.size();i++){
		fs<<control::site_name[i]<<"\t";
	}
	fs<<std::endl;
	for(size_t i=0;i<control::site_name.size();i++){
		fs<<std::setw(8)<<std::fixed<<control::vm[tick]<<"\t";
		tick++;
	}
	fs<<"(the final site parameter change range is zero due to Charge-Neutral)";
	fs<<std::endl;
}
void write_defile(int databasetick){
	std::fstream fs;
	fs.open(control::deopt[databasetick].c_str(),std::fstream::out);
	fs<<"#Difference(ev)"<<"\t"<<"MD(ev)"<<"\t"<<"DFT(ev)"<<std::endl;
  int size=0;
	size=control::ionsize[databasetick];
	for(size_t i=0;i<size;i++){
		fs<<std::setw(8)<<std::fixed<<control::diffenergy[databasetick][i]<<"\t"<<control::mdenergy[databasetick][i]<<"\t"<<control::dftenergy[databasetick][i]<<std::endl;
	}
	fs.close();
}
