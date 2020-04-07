#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <sstream>
#include <map>
#include <cstdlib>
#include <cmath>
std::vector<std::vector<double> > extract_periodical(std::string dftinput){
	std::fstream fs;
	fs.open(dftinput.c_str(),std::fstream::in);
	std::vector<std::vector<double> > periodical(3,std::vector<double>(3,0.0));
	std::string temp;
	std::stringstream temp_stream;
	while(getline(fs,temp)){
		if(temp.find("CELL_PARAMETERS")!=std::string::npos){
			for(int i=0;i<3;i++){
				getline(fs,temp);
				temp_stream.clear();
				temp_stream.str(temp);
				for(int j=0;j<3;j++){
					temp_stream>>periodical[i][j];
				}
				temp_stream.clear();
			}
		}
	}
	fs.close();
	return periodical;
}
double vectordot(std::vector<double>& a,std::vector<double>& b,int c){
	double sum=0.0;
	for(int i=0;i<c;i++){
		sum=sum+a[i]*b[i];
	}
	return sum;
}
void normalize(std::vector<double>& input,int c){
	double sum=0.0;
	for(int i=0;i<3;i++){
		sum=sum+input[i]*input[i];
	}
	for(int i=0;i<3;i++){
		input[i]=input[i]/sqrt(sum);
	}
}
std::vector<std::vector<double> > changeaxis(std::vector<std::vector<double> > &input){
	std::vector<double> x(3,0.0);
	std::vector<double> y(3,0.0);
	std::vector<double> z(3,0.0);
	std::vector<std::vector<double> > axis(3,std::vector<double>(3,0.0));
	axis[0]=input[0];
	normalize(axis[0],3);
	for(int i=0;i<3;i++){
		axis[1][i]=input[1][i]-vectordot(axis[0],input[1],3)*axis[0][i];
	}
	normalize(axis[1],3);
	for(int i=0;i<3;i++){
	  axis[2][i]=input[2][i]-vectordot(axis[0],input[2],3)*axis[0][i]-vectordot(axis[1],input[2],3)*axis[1][i];
	}
	normalize(axis[2],3);
	std::vector<std::vector<double> > re(3,std::vector<double>(3,0.0));
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
		 re[i][j]=vectordot(input[i],axis[j],3);
		}
	}
	return re;
}
std::vector<double> coordinate(std::vector<double>& fraction,std::vector<std::vector<double> >& axis){
	std::vector<double> re(3,0.0);
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			re[j]=re[j]+axis[i][j]*fraction[i];
		}
	}
	return re;
}
std::vector<double> coordinate(double nx,double ny,double nz,std::vector<std::vector<double> >& axis){
	std::vector<double> re(3,0.0);
	std::vector<double> fraction(3,0.0);
	fraction[0]=nx;
	fraction[1]=ny;
	fraction[2]=nz;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			re[j]=re[j]+axis[i][j]*fraction[i];
		}
	}
	return re;
}
int main(){
	std::vector<std::vector<double> > period=extract_periodical("bto.in");
	std::vector<std::vector<double> > axis=changeaxis(period);
	std::vector<double> coord(3,0.0);
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			axis[i][j]=axis[i][j]/2;
	double asitecharge=1.446816;
	double bsitecharge=1.34129;
	double ositecharge=-0.929368;
	int cell=10;
	int cellindft=2;
	std::fstream fs;
	fs.open("mixdata.BTO",std::fstream::out);
	fs<<"#LAMMPS 100 water\n";
	fs<<"\n";
	fs<<cell*cell*cell*5<<" atoms\n";
	fs<<"3 atom types\n";
	fs<<"\n";
	fs<<"0.0 "<<std::setprecision(11)<<axis[0][0]*cell<<" xlo xhi\n";
	fs<<"0.0 "<<std::setprecision(11)<<axis[1][1]*cell<<" ylo yhi\n";
	fs<<"0.0 "<<std::setprecision(11)<<axis[2][2]*cell<<" zlo zhi\n";
	fs<<axis[1][0]*cell<<" "<<axis[2][0]*cell<<" "<<axis[2][1]*cell<<" xy xz yz\n";
	fs<<"\n";
	fs<<"\n";
	fs<<"Masses\n";
	fs<<"\n";
	fs<<"1 208.98 #Bi\n";
	fs<<"2 55.845 #Fe\n";
	fs<<"3 15.9999 #O\n";
	fs<<"\n";
	fs<<"Atoms\n";
	fs<<"\n";
	for(int k=0;k<cell;k++)
		for(int j=0;j<cell;j++)
			for(int i=0;i<cell;i++){
			coord=coordinate(i,j,k,axis);
			fs<<i+j*cell+k*cell*cell+1<<" "<<1<<" "<<1<<" "<<asitecharge<<" "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<std::endl;
			}
	for(int k=0;k<cell;k++)
		for(int j=0;j<cell;j++)
			for(int i=0;i<cell;i++){
			coord=coordinate(i+0.250196608,j+0.277942506,k+0.285637675,axis);
			fs<<i+j*cell+k*cell*cell+1+cell*cell*cell<<" "<<1<<" "<<2<<" "<<bsitecharge<<" "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<std::endl;
			}
	for(int k=0;k<cell;k++)
		for(int j=0;j<cell;j++)
			for(int i=0;i<cell;i++){
			coord=coordinate(i+0.25,j+0.295163605,k+0.088030717,axis);
			fs<<i+j*cell+k*cell*cell+1+cell*cell*cell*2<<" "<<1<<" "<<3<<" "<<ositecharge<<" "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<std::endl;
			}
	for(int k=0;k<cell;k++)
		for(int j=0;j<cell;j++)
			for(int i=0;i<cell;i++){
			coord=coordinate(i+0.25,j+0.027278961,k+0.351475630,axis);
			fs<<i+j*cell+k*cell*cell+1+cell*cell*cell*3<<" "<<1<<" "<<3<<" "<<ositecharge<<" "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<std::endl;
			}
	for(int k=0;k<cell;k++)
		for(int j=0;j<cell;j++)
			for(int i=0;i<cell;i++){
			coord=coordinate(i+0.000000221,j+0.278355902,k+0.355330868,axis);
			fs<<i+j*cell+k*cell*cell+1+cell*cell*cell*4<<" "<<1<<" "<<3<<" "<<ositecharge<<" "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<std::endl;
			}
}
