#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sstream>
#include <map>
#include <cstdlib>
/*this application use to use the opt.out file to run MD simulations and extract the MD simulation result to make another iteration dft input*/
void generatelammpsinput(std::string opt,int atom,std::map<std::string,double>& chargemap){
	std::stringstream temp_stream;
	std::fstream fs;
	fs.open(opt.c_str(),std::fstream::in);
	std::string temp;
	std::vector<std::string> species(atom,"");
	std::vector<double> charge(atom,0.0);
	std::vector<std::vector<double> > parameter(atom*(atom+1)/2,std::vector<double>(12,0.00));
	while(getline(fs,temp)){
		if(temp.find("new-optimized")!=std::string::npos){
			chargemap.clear();
			for(size_t i=0;i<atom*(atom+1)/2;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				for(size_t j=0;j<12;j++){
					temp_stream>>parameter[i][j];
				}
				temp_stream.clear();
			}
			getline(fs,temp);
			temp_stream.clear();
			temp_stream.str(temp);
			for(size_t i=0;i<atom;i++){
				temp_stream>>species[i];
			}
			getline(fs,temp);
			temp_stream.clear();
			temp_stream.str(temp);
			for(size_t i=0;i<atom;i++){
				temp_stream>>charge[i];
			}
			for(size_t i=0;i<atom;i++){
			chargemap.insert(std::pair<std::string,double>(species[i],charge[i]));
			}
		}
	}
	fs.close();
	fs.open("input.lammps",std::fstream::out);
	fs<<"# temperarory input for lammps"<<std::endl;
	fs<<" "<<std::endl;
	fs<<"units     metal"<<std::endl;
	fs<<"atom_style full"<<std::endl;
	fs<<"boundary p p p"<<std::endl;
	fs<<"kspace_style pppm 1.0e-4"<<std::endl;
	fs<<"pair_style hybrid/overlay  12lj/cut/coul/long 8.0 8.0 bv 2.0 8.0 bvv 2.0 8.0"<<std::endl;
	fs<<"angle_style harmonic"<<std::endl;
	fs<<"read_data mixdata.BTO"<<std::endl;
	fs<<"#read_restart BTO.restart"<<std::endl;
	for(size_t i=0;i<3;i++){
		fs<<std::endl;
	}
	size_t tick=0;
	for(size_t i=1;i<=atom;i++){
		for(size_t j=i;j<=atom;j++){
			fs<<"pair_coeff"<<" "<<i<<" "<<j<<" "<<"12lj/cut/coul/long 2.0"<<" "<<std::fixed<<std::setprecision(7)<<parameter[tick][7]<<std::endl;
			tick++;
		}
	}
	fs<<std::endl;
	fs<<std::endl;
	fs<<"#                   r0   Nij    S     v0 rcut"<<std::endl;
	tick=0;
	for(size_t i=1;i<=atom;i++){
		for(size_t j=i;j<=atom;j++){
			fs<<"pair_coeff"<<" "<<i<<" "<<j<<" "<<"bv"<<" "<<std::fixed<<std::setprecision(7)<<parameter[tick][0]<<" "<<parameter[tick][1]<<" "<<parameter[tick][2]<<" "<<parameter[tick][4]<<" "<<parameter[tick][5]<<std::endl;
			tick++;
		}
	}
	fs<<"#                    r0  Nij     S     Bvv0  rcut"<<std::endl;
	tick=0;
	for(size_t i=1;i<=atom;i++){
		for(size_t j=i;j<=atom;j++){
			fs<<"pair_coeff"<<" "<<i<<" "<<j<<" "<<"bvv"<<" "<<std::fixed<<std::setprecision(7)<<parameter[tick][0]<<" "<<parameter[tick][1]<<" "<<parameter[tick][10]<<" "<<parameter[tick][11]<<" "<<parameter[tick][5]<<std::endl;
			tick++;
		}
	}
	fs<<std::endl;
	fs<<std::endl;
	fs<<"neighbor        2.0 bin"<<std::endl;
	fs<<"neigh_modify    every 1"<<std::endl;
	fs<<"# time unit ps"<<std::endl;
	fs<<"timestep         0.001"<<std::endl;
	fs<<std::endl;
	fs<<"#group Pb id 1:512"<<std::endl;
	fs<<"#group Ti id 513:1024"<<std::endl;
	fs<<"#group O1 id 1025:1536"<<std::endl;
	fs<<"#group O2 id 1537:2048"<<std::endl;
	fs<<"#group O3 id 2049:2560"<<std::endl;
	fs<<"thermo          100"<<std::endl;
	fs<<"thermo_style custom step temp eangle etotal press vol lx ly lz"<<std::endl;
	fs<<"thermo_modify line one format float %12.5f"<<std::endl;
	fs<<std::endl;
	fs<<std::endl;
	fs<<"minimize 1.0e-6 1.0e-8 1000 10000"<<std::endl;
	fs<<std::endl;
	fs<<"fix             2 all npt temp 10.0 10.0 1.0 aniso 1.01325 1.01325 5.0"<<std::endl;
	fs<<"run 60000"<<std::endl;
	fs<<"unfix 2"<<std::endl;
	fs<<std::endl;
	fs<<"fix 4 all npt temp 10.0 10.0 1.0 aniso 1.01325 1.01325 5.0"<<std::endl;
	fs<<"dump 4 all custom 2000 dump.xyz x y z"<<std::endl;
	fs<<"dump_modify 4 sort id"<<std::endl;
	fs<<"run 10000"<<std::endl;
	fs<<"unfix 4"<<std::endl;
	fs<<std::endl;
	fs.close();
}
std::vector<std::string> extractspecies(std::string dftinput,int atomnum){
	std::fstream fs;
	fs.open(dftinput.c_str(),std::fstream::in);
	std::vector<std::string> species;
	std::string temp;
	std::stringstream temp_stream;
	std::string single_species;
	while(getline(fs,temp)){
		if(temp.find("ATOMIC_POSITIONS")!=std::string::npos){
			for(int i=0;i<atomnum;i++){
			getline(fs,temp);
			temp_stream.clear();
			temp_stream.str(temp);
			temp_stream>>single_species;
			species.push_back(single_species);
			}
		}
	}
	fs.close();
	return species;
}
std::vector<std::vector<double> > extract_starting_position(std::string dftinput,int atomnum){
	std::fstream fs;
	fs.open(dftinput.c_str(),std::fstream::in);
	std::vector<double> position(3,0.0);
	std::vector<double> periodical(3,0.0);
	std::vector<std::vector<double> > allposition;
	std::string species;
	std::string temp;
	std::stringstream temp_stream;
	while(getline(fs,temp)){
		if(temp.find("CELL_PARAMETERS")!=std::string::npos){
			for(int i=0;i<3;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				for(int j=0;j<=i;j++){
					temp_stream>>periodical[i];
				}
				temp_stream.clear();
			}
		}
		if(temp.find("ATOMIC_POSITIONS")!=std::string::npos){
			if(temp.find("angstrom")!=std::string::npos){
				for(int i=0;i<atomnum;i++){
					getline(fs,temp);
					temp_stream.clear();
					temp_stream.str(temp);
					temp_stream>>species;
					for(size_t j=0;j<3;j++){
						temp_stream>>position[j];
					}
					temp_stream.clear();
					allposition.push_back(position);
					}
				}
			else if(temp.find("crystal")!=std::string::npos){
				for(int i=0;i<atomnum;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				temp_stream>>species;
				for(size_t j=0;j<3;j++){
					temp_stream>>position[j];
					position[j]=position[j]*periodical[j];
				}
				temp_stream.clear();
				allposition.push_back(position);
			}		   
			}
		}
	}
	fs.close();
	return allposition;
}
std::vector<double> extract_periodical(std::string dftinput,int atomnum){
	std::fstream fs;
	fs.open(dftinput.c_str(),std::fstream::in);
	std::vector<double> periodical(3,0.0);
	std::string species;
	std::string temp;
	std::stringstream temp_stream;
	while(getline(fs,temp)){
		if(temp.find("CELL_PARAMETERS")!=std::string::npos){
			for(int i=0;i<3;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				for(int j=0;j<=i;j++){
					temp_stream>>periodical[i];
				}
				temp_stream.clear();
			}
		}
	}
	fs.close();
	return periodical;
}
void generatelammpsdata(std::string dftinput,int atomnum,std::map<std::string,int>& speciesmap, std::map<std::string,double>& massmap,std::map<std::string,double>& chargemap){
	std::fstream fs;
	std::vector<double> periodical=extract_periodical(dftinput,atomnum);
	std::vector<std::vector<double> > allposition=extract_starting_position(dftinput,atomnum);
	std::vector<std::string> species=extractspecies(dftinput,atomnum);
	fs.open("mixdata.BTO",std::fstream::out);
	fs<<"#LAMMPS 100 water\n";
	fs<<"\n";
	fs<<atomnum<<" atoms\n";
	fs<<"4 atom types\n";
  fs<<"\n";
	fs<<"0.0 "<<std::setprecision(11)<<periodical[0]<<" "<<"xlo xhi\n";
	fs<<"0.0 "<<std::setprecision(11)<<periodical[1]<<" "<<"ylo yhi\n";
	fs<<"0.0 "<<std::setprecision(11)<<periodical[2]<<" "<<"zlo zhi\n";
	fs<<"\n";
	fs<<"\n";
	fs<<"Masses\n";
	fs<<"\n";
	for(std::map<std::string,int>::iterator a=speciesmap.begin();a!=speciesmap.end();a++){
		fs<<a->second<<" "<<massmap[a->first]<<" #"<<a->first<<"\n";
	}
	fs<<"\n";
	fs<<"\n";
	fs<<"Atoms\n";
	fs<<"\n";
	for(int i=0;i<atomnum;i++){
		fs<<i+1<<" "<<1<<" "<<speciesmap[species[i]]<<" "<<std::setprecision(9)<<chargemap[species[i]]<<" ";
		for(size_t j=0;j<3;j++){
			fs<<std::setprecision(11)<<allposition[i][j]<<" ";
		}
		fs<<"\n";
	}
	fs.close();
}
void generatedft(std::string dftinput,int atomnum){
	std::vector<std::string> species=extractspecies(dftinput,atomnum);
	std::vector<std::string> head;
	std::vector<std::string> dumpfile;
	std::fstream fs;
	fs.open(dftinput.c_str(),std::fstream::in);
	std::string temp;
	double periodical[3]={0.0,0.0,0.0};
	std::stringstream temp_stream;
	double x1,x2;
	while(getline(fs,temp)){
		if(temp.find("CELL_PARAMETERS")!=std::string::npos){
		break;
		}
		else{
			head.push_back(temp);
		}
	}
	fs.close();
	fs.open("dump.xyz",std::fstream::in);
	while(getline(fs,temp)){
		dumpfile.push_back(temp);
	}
	fs.close();
	std::fstream fsdftout;
	std::string dftoutfile;
	std::stringstream ss;
	for(size_t i=0;i<5;i++){
		ss.str("");
		ss<<i;
		dftoutfile="bto.in"+ss.str();
		fsdftout.open(dftoutfile.c_str(),std::fstream::out);
		for(std::vector<std::string>::iterator a=head.begin();a!=head.end();a++){
			fsdftout<<*a<<std::endl;
		}
		for(size_t j=0;j<3;j++){
			temp=dumpfile[i*(atomnum+9)+5+j];
			temp_stream.clear();
			temp_stream.str(temp);
			temp_stream>>x1;
			temp_stream>>x2;
			periodical[j]=x2-x1;
		}
		fsdftout<<"CELL_PARAMETERS (angstrom)"<<std::endl;
		for(size_t j=0;j<3;j++){
			for(size_t k=0;k<3;k++){
				if(j==k){
					fsdftout<<periodical[j]<<" ";
				}
				else{
					fsdftout<<0.0<<" ";
				}
			}
			fsdftout<<std::endl;
		}
		fsdftout<<"ATOMIC_POSITIONS (angstrom)"<<std::endl;
		for(size_t j=0;j<atomnum;j++){
			fsdftout<<species[j]<<" "<<dumpfile[i*(atomnum+9)+9+j]<<std::endl;
		}
		fsdftout.close();
	}
}
int main(){
	/*generate lammps input*/
	int atomtype=4;
	int atomnum=40;
	std::map<std::string,int> speciesmap;
	std::map<std::string,double> massmap;
	speciesmap.insert(std::pair<std::string,int>("Ca",2));
	speciesmap.insert(std::pair<std::string,int>("Ba",1));
	speciesmap.insert(std::pair<std::string,int>("Zr",3));
	speciesmap.insert(std::pair<std::string,int>("O",4));
	massmap.insert(std::pair<std::string,double>("Ca",40.078));
	massmap.insert(std::pair<std::string,double>("Ba",137.33));
	massmap.insert(std::pair<std::string,double>("Zr",91.224));
	massmap.insert(std::pair<std::string,double>("O",15.999));
	std::map<std::string,double> chargemap;
	/*generate lammps input file*/
	generatelammpsinput("opt.out",atomtype,chargemap);
	/*generate lammps datafile*/
  generatelammpsdata("bto.in",atomnum,speciesmap,massmap,chargemap);
	/*run the lammps file*/
	std::string lammps_command="./lmp_serial < input.lammps";
	std::system(lammps_command.c_str());
	/*collect the trajectories to dft input file*/
	generatedft("bto.in",atomnum);
}
