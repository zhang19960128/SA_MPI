#include <fstream>
#include <string>
#include "atom.h"
#include <sstream>
#include <iostream>
#include <stdlib.h>
double** getparameter(std::string pair_style,std::string file){
      std::fstream fs;
      fs.open(file,std::fstream::in);
      std::string line;
      std::istringstream stream1;
      std::string temp;
      int count=0;
      do{
	    getline(fs,line);
	    stream1.str(line);
	    while(stream1>>temp){
               if(temp=="pair_coeff"){
                  stream1>>temp;
                  stream1>>temp;
                  stream1>>temp;
                  if(temp==pair_style){
                     count++;
                  }
               }
         }
         stream1.clear();
	}while(!fs.eof());
        /*test how many atoms are in the input*/
        size_t type=1;
        for(type=1;type<30;type++){
            if(type*(type+1)/2==count){
               break;
            }
        }
        fs.close();
        if(type==30){
           std::cout<<"error, input pair loss or exceed for pair potential: "<<pair_style<<std::endl;
           exit(EXIT_FAILURE);
        }
        else{
          std::cout<<"In this ### "<<pair_style<<" ###: pair potential, we found " <<type<<" types of atoms"<<std::endl;
        }
      double** para=new double* [count];
      fs.open(file,std::fstream::in);
      int tick1,tick2;
      int onedim;
      int* base=new int[type];
      for(size_t j=0;j<type;j++){
         base[j]=(type-j+type)*(j+1)/2;
      }
      do{
	    getline(fs,line);
	    stream1.str(line);
	    while(stream1>>temp){
               if(temp=="pair_coeff"){
                  stream1>>tick1;
                  stream1>>tick2;
                  if(tick1>tick2){
                     std::cout<<"when input the pair parameter, please make the second type tick bigger than first"<<std::endl;
                     exit(EXIT_FAILURE);
                  }
                  stream1>>temp;
                  onedim=(tick1-1!=0)*base[(tick1-2)>0?(tick1-2):0]+(tick2-tick1);
                  if(temp==pair_style){
                  if(temp=="12lj/cut/coul/long"){
                     para[onedim]=new double[2];
                     stream1>>para[onedim][0];
                     stream1>>para[onedim][1];
                  }
                  else if(temp=="bv"){
                     para[onedim]=new double[5];
                     for(size_t i=0;i<5;i++){
                        stream1>>para[onedim][i];
                     }
                  }
                  else if(temp=="bvv"){
                     para[onedim]=new double[5];
                     for(size_t i=0;i<5;i++){
                        stream1>>para[onedim][i];
                     }
                  }
                  else{
                     std::cout<<"unknow pair style"<<std::endl;
                     exit(EXIT_FAILURE);
                  }
                }
               }
         }
         stream1.clear();
	}while(!fs.eof());
        fs.close();
	return para;
}
atom* configuration(std::string file){
	std::fstream fs;
	fs.open(file,std::fstream::in);
	std::string line;
	std::string one;
	std::string two;
	std::istringstream stream1;
	int atomnum;
	int atomtype;
	do{
		getline(fs,line);
	  stream1.str(line);
		two="";
		while(stream1>>one){
			if(one=="atoms"&&two.length()!=0){
				atomnum=std::stoi(two);
			}
			two=one;
		}
		stream1.clear();
	}while(!fs.eof());
	fs.close();
	/*we now know how many atoms are in the configuration file*/
	atom* config=new atom[atomnum];
	fs.open(file,std::fstream::in);
	int ticknum,groupID,type;
	double charge,x,y,z;
	do{
		getline(fs,line);
		if(line=="Atoms"){
			getline(fs,line);
			while(line.length()==0){
				getline(fs,line);
			}
			for(size_t j=0;j<atomnum;j++){
				stream1.str(line);
				stream1>>ticknum;
				stream1>>groupID;
				stream1>>type;
				stream1>>charge;
				stream1>>x;
				stream1>>y;
				stream1>>z;
				config[ticknum-1].charge=charge;
				config[ticknum-1].type=type-1;
				config[ticknum-1].position[0]=x;
				config[ticknum-1].position[1]=y;
				config[ticknum-1].position[2]=z;
				stream1.clear();
				getline(fs,line);
			}
		}
	}while(!fs.eof());
	return config;
}
