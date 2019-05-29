#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
double bohr=0.529177249;
double ha=27.2113845;
double Ry=13.60569253;
double gpa=0.1;
void write(std::fstream& out,int atomnum,std::vector<std::vector<double> >& position,std::vector<std::vector<double> >& v3,std::vector<std::vector<double> >& force,double& total_energy,std::vector<std::vector<double> >& pressure,double weight,int& tick){
		out<<"******* Reduced ionic position : "<<tick<<std::endl;
		for(size_t i=0;i<atomnum;i++){
			for(size_t j=0;j<3;j++){
				out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<position[i][j]<<"\t";
			}
			 out<<std::endl;
		}
		out<<"******* Lattice unit vectors"<<std::endl;
		for(size_t i=0;i<3;i++){
			out<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][0]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][1]<<"\t"<<std::fixed<<std::setprecision(15)<<std::setw(15)<<v3[i][2]<<std::endl;
		}
		out<<"******* Lattice lengths (A)"<<std::endl;
		out<<"1.000 1.000 1.000"<<std::endl;
		out<<"******* Ionic forces (eV/A)"<<std::endl;
		for(size_t i=0;i<atomnum;i++){
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
		out<<weight<<std::endl;
		std::cout<<tick++<<std::endl;

}
int main(int argc,char* argv[]){
	int interval=1;
	std::string dftoutfile=argv[1];
	std::string outfile=argv[2];
	int natom=40;
	/*units that are useful*/
	std::fstream dftout;
	dftout.open(dftoutfile.c_str(),std::fstream::in);
	std::fstream out;
	out.open(outfile.c_str(),std::fstream::out);
	std::string temp;
	int start=std::atoi(argv[3]);
	int end=start+8000;
	std::istringstream stream1;
	std::string substr;
	double posit;
	std::vector<double> v1(3,0.0);
	std::vector<std::vector<double> > v3(3,v1);
	std::vector<std::vector<double> > pressure(3,v1);
	std::vector<std::vector<double> > force(natom,v1);
	std::vector<std::vector<double> > position(natom,v1);
	double weight=1.0;
  for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++){
          if(i!=j)
              v3[i][j]=0;
      }
  v3[0][0]=8.9258090014;
  v3[1][1]=8.9258090014;
  v3[2][2]=12.623000144;
	double total_energy;
	int write_signal=1;
	int P_signal=0;
	int F_signal=0;
	int E_signal=0;
	int S_signal=0;
	do{
		getline(dftout,temp);
		P_signal=0;
		F_signal=0;
		E_signal=0;
		S_signal=0;
		if(temp.find("CELL_PARAMETERS (angstrom)")!=std::string::npos){
			for(size_t i=0;i<3;i++){
			getline(dftout,temp);
			stream1.str(temp);
			for(size_t j=0;j<3;j++){
				stream1>>v3[i][j];
			}
			stream1.clear();
			}
		}
		else if(temp.find("ATOMIC_POSITIONS (angstrom)")!=std::string::npos){
			P_signal=1;
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
				position[i][j]=posit;
			 }
			 stream1.clear();
			}
		do{
						getline(dftout,temp);
		}while(temp.find("!    total energy")==std::string::npos && !dftout.eof());
		if(temp.find("!    total energy")!=std::string::npos){
		stream1.clear();
	  stream1.str(temp);
	  stream1>>substr;//!
	  stream1>>substr;//total
	  stream1>>substr;//energy
	  stream1>>substr;//=
	  stream1>>total_energy;//total energy
	  stream1.clear();
		E_signal=1;
		}
		do{
						getline(dftout,temp);
		}while(temp.find("Forces acting on atoms")==std::string::npos && !dftout.eof());
		if(temp.find("Forces acting on atoms")!=std::string::npos){	
			F_signal=1;
			getline(dftout,temp);
			for(size_t i=0;i<natom;i++){
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
		do{
						getline(dftout,temp);
		}while(temp.find("total   stress")==std::string::npos && !dftout.eof());
		if(temp.find("total   stress")!=std::string::npos){
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
			S_signal=1;
		}
		if(F_signal*E_signal*S_signal*P_signal==1){
        write(out,natom,position,v3,force,total_energy,pressure,weight,start);
		}
			if(start==end) break;
		}
		else if(temp=="ATOMIC_POSITIONS (crystal)"){
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
				position[i][j]=posit;
			 }
			 stream1.clear();
		}
		P_signal=1;
		do{
						getline(dftout,temp);
		}while(temp.find("!    total energy")==std::string::npos && !dftout.eof());
		if(temp.find("!    total energy")!=std::string::npos){
		stream1.clear();
	  stream1.str(temp);
	  stream1>>substr;//!
	  stream1>>substr;//total
	  stream1>>substr;//energy
	  stream1>>substr;//=
	  stream1>>total_energy;//total energy
	  stream1.clear();
		E_signal=1;
		}
		do{
						getline(dftout,temp);
		}while(temp.find("Forces acting on atoms")==std::string::npos && !dftout.eof());
		if(temp.find("Forces acting on atoms")!=std::string::npos){	
		        getline(dftout,temp);
			for(size_t i=0;i<natom;i++){
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
			F_signal=1;
		}
		do{
						getline(dftout,temp);
		}while(temp.find("total   stress")==std::string::npos && !dftout.eof());
		if(temp.find("total   stress")!=std::string::npos){
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
		S_signal=1;
		}
		if(F_signal*E_signal*S_signal*P_signal==1){
        write(out,natom,position,v3,force,total_energy,pressure,weight,start);
		}
		}
		if(start==end) break;
	}while(!dftout.eof());
	out.close();
}
