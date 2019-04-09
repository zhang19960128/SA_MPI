#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
int main(int argc,char* argv[]){
	int interval=1;
	std::string dftoutfile=argv[1];
	std::string outfile="IonFor.dat.MD.jiahao";
	int natom=40;
	int nline=41;
	/*units that are useful*/
	double bohr=0.529177249;
	double ha=27.2113845;
	double Ry=13.60569253;
	double gpa=0.1;
	std::fstream dftout;
	dftout.open(dftoutfile.c_str(),std::fstream::in);
	std::fstream out;
	out.open(outfile.c_str(),std::fstream::out);
	std::string temp;
	int start=0;
	int end=start+3000;
	std::istringstream stream1;
	std::string substr;
	double posit;
	std::vector<double> v1(3,0.0);
	std::vector<std::vector<double> > v3(3,v1);
	std::vector<std::vector<double> > pressure(3,v1);
	std::vector<std::vector<double> > force(natom,v1);
	double total_energy;
	do{
		getline(dftout,temp);
		if(temp.find("!    total energy")!=std::string::npos){
			stream1.str(temp);
			stream1>>substr;//!
			stream1>>substr;//total
			stream1>>substr;//energy
			stream1>>substr;//=
			stream1>>total_energy;//total energy
		  stream1.clear();
		}
		else if(temp.find("Forces acting on atoms (cartesian axes, Ry/au):")!=std::string::npos){
			getline(dftout,temp);
			for(size_t i=0;i<40;i++){
			getline(dftout,temp);
			stream1.str(temp);
			stream1>>substr;//atom
			stream1>>substr;//number
			stream1>>substr;//type
			stream1>>substr;//number
			stream1>>substr;//force
			stream1>>substr;//=
		  for(size_t j=0;j<3;j++){
				stream1>>force[i][j];
			}
			stream1.clear();
			}
				}
		else if(temp.find("total   stress") != std::string::npos){
			for(size_t i=0;i<3;i++){
				getline(dftout,temp);
				stream1.str(temp);
				for(size_t j=0;j<3;j++){
					stream1>>posit;//get ride of the non-useful stress double.
				}
				for(size_t j=0;j<3;j++){
					stream1>>pressure[i][j];
				}
				stream1.clear();
			}
		}
		else if(temp=="CELL_PARAMETERS (angstrom)"){
			for(size_t i=0;i<3;i++){
			getline(dftout,temp);
			stream1.str(temp);
			for(size_t j=0;j<3;j++){
				stream1>>v3[i][j];
			}
			stream1.clear();
			}
		}
		else if(temp=="ATOMIC_POSITIONS (crystal)"){
		out<<"******* Reduced ionic position : "<<start<<std::endl;
		for(size_t i=0;i<natom;i++){
			 getline(dftout,temp);
			 stream1.str(temp);
			 stream1>>substr;
			 for(size_t j=0;j<3;j++){
			 	stream1>>posit;
				if(posit<0){
					posit=posit+1;
				}
				else if(posit>1){
					posit=posit-1;
				}
				else posit=posit;
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<posit<<"\t";
			 }
			 out<<std::endl;
			 stream1.clear();
		}
out<<"******* Lattice unit vectors"<<std::endl;
		for(size_t i=0;i<3;i++){
		out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][0]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][1]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][2]<<std::endl;
		}
		out<<"******* Lattice lengths (A)"<<std::endl;
		out<<"1.000 1.000 1.000"<<std::endl;
		out<<"******* Ionic forces (eV/A)"<<std::endl;
		for(size_t i=0;i<natom;i++){
			for(size_t j=0;j<3;j++){
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<force[i][j]*Ry/bohr<<"\t";
			}
			out<<std::endl;
		}
		out<<"******* Total energy (eV/supercell)"<<std::endl;
    out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<total_energy*Ry<<std::endl;
		out<<"******* Stress (GPa)"<<std::endl;
		for(size_t i=0;i<3;i++){
			for(size_t j=0;j<3;j++){
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<pressure[i][j]*gpa<<"\t";
			}
		    out<<std::endl;
		}
		out<<"******* Weight "<<std::endl;
		out<<"1.0"<<std::endl;
		std::cout<<start++<<std::endl;
		if(start==end) break;
		}
		else if(temp=="ATOMIC_POSITIONS (angstrom)"){
			out<<"******* Reduced ionic position : "<<start<<std::endl;
	  	for(size_t i=0;i<natom;i++){
			 getline(dftout,temp);
			 stream1.str(temp);
			 stream1>>substr;
			 for(size_t j=0;j<3;j++){
			 	stream1>>posit;
				posit=posit/v3[j][j];
				if(posit<0){
					posit=posit+1;
				}
				else if(posit>1){
					posit=posit-1;
				}
				else posit=posit;
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<posit<<"\t";
			 }
			 out<<std::endl;
			 stream1.clear();
			}
		out<<"******* Lattice unit vectors"<<std::endl;
		for(size_t i=0;i<3;i++){
		out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][0]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][1]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][2]<<std::endl;
		}
		out<<"******* Lattice lengths (A)"<<std::endl;
		out<<"1.000 1.000 1.000"<<std::endl;
		out<<"******* Ionic forces (eV/A)"<<std::endl;
		for(size_t i=0;i<natom;i++){
			for(size_t j=0;j<3;j++){
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<force[i][j]*Ry/bohr<<"\t";
			}
			out<<std::endl;
		}
		out<<"******* Total energy (eV/supercell)"<<std::endl;
    out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<total_energy*Ry<<std::endl;
		out<<"******* Stress (GPa)"<<std::endl;
		for(size_t i=0;i<3;i++){
			for(size_t j=0;j<3;j++){
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<pressure[i][j]*gpa<<"\t";
			}
		    out<<std::endl;
		}
		out<<"******* Weight "<<std::endl;
		out<<"1.0"<<std::endl;
		std::cout<<start++<<std::endl;
		if(start==end) break;
		}
	}while(!dftout.eof());
	out.close();
}
