#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sstream>
int main(int argc,char* argv[]){
	std::stringstream temp_stream;
	int atom;
	temp_stream.str(argv[1]);
	temp_stream>>atom;
	std::fstream fs;
	fs.open("opt.out",std::fstream::in);
	std::string temp;
	std::vector<std::vector<double> > parameter(atom*(atom+1)/2,std::vector<double>(12,0.00));
	while(getline(fs,temp)){
		if(temp.find("new-optimized")!=std::string::npos){
			for(size_t i=0;i<atom*(atom+1)/2;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				for(size_t j=0;j<12;j++){
					temp_stream>>parameter[i][j];
				}
				temp_stream.clear();
			}
		}
	}
	fs.close();
	for(size_t i=0;i<atom*(atom+1)/2;i++){
		for(size_t j=0;j<12;j++){
			std::cout<<parameter[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
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
	fs<<"fix 1 all nvt temp 100.0 100.0 1.0"<<std::endl;
	fs<<"run 50000"<<std::endl;
	fs<<"unfix 1"<<std::endl;
	fs<<std::endl;
	fs<<"fix 2 all npt temp 100.0 100.0 1.0 aniso 1.01325 1.01325 5.0"<<std::endl;
	fs<<"run 50000"<<std::endl;
	fs<<"unfix 2"<<std::endl;
	fs<<std::endl;
		fs<<"fix 3 all nvt temp 200.0 200.0 1.0"<<std::endl;
	fs<<"run 50000"<<std::endl;
	fs<<"unfix 3"<<std::endl;
	fs<<std::endl;
	fs<<"fix 4 all npt temp 200.0 200.0 1.0 aniso 1.01325 1.01325 5.0"<<std::endl;
	fs<<"dump 4 all custom 200 dump.xyz x y z"<<std::endl;
	fs<<"dump_modify 4 sort id"<<std::endl;
	fs<<"run 250000"<<std::endl;
	fs<<"unfix 4"<<std::endl;
	fs<<std::endl;
}
