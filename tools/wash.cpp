#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
int main(){
	std::vector<std::vector<double> > period(3,std::vector<double>(2,0.0));
	std::string boxflag="ITEM: BOX BOUNDS pp pp pp";
	std::string positflag="ITEM: ATOMS x y z";
	std::vector<std::vector<double> > before(40,std::vector<double>(3,0.0));
	std::vector<std::vector<double> > now(40,std::vector<double>(3,0.0));
	std::fstream fs,fsout;
	fs.open("dump.xyz.E",std::fstream::in);
	fsout.open("traject",std::fstream::out);
	std::string temp;
	std::istringstream stream1;
	stream1.clear();
	int atom=40;
	int frame=0;
	while(getline(fs,temp)){
		if(temp.find(boxflag)!=std::string::npos){
			for(size_t i=0;i<3;i++){	
				getline(fs,temp);
				stream1.clear();
				stream1.str(temp);
				stream1>>period[i][0];
				stream1>>period[i][1];
			}
		}
		if(temp.find(positflag)!=std::string::npos){
			for(size_t i=0;i<atom;i++){
				for(size_t j=0;j<3;j++){
					before[i][j]=now[i][j];
				}
			}
			for(size_t i=0;i<atom;i++){
				getline(fs,temp);
				stream1.clear();
				stream1.str(temp);
				for(size_t j=0;j<3;j++){
					stream1>>now[i][j];
				}
			}
		if(frame==0){
			for(size_t i=0;i<atom;i++){
				for(size_t j=0;j<3;j++){
					before[i][j]=now[i][j];
				}
			}
		}
		frame++;
		/*starting to wash the data*/
		for(size_t i=0;i<atom;i++){
			for(size_t j=0;j<3;j++){
				std::cout<<now[i][j]-before[i][j]<<std::endl;
				if(now[i][j]-before[i][j]>(period[j][1]-period[j][0])/2.0){
					now[i][j]=now[i][j]-(period[j][1]-period[j][0]);
					std::cout<<i<<" "<<j<<std::endl;
				}
				else if(now[i][j]-before[i][j]<-1*(period[j][1]-period[j][0])/2.0){
					now[i][j]=now[i][j]+(period[j][1]-period[j][0]);
					std::cout<<i<<" "<<j<<std::endl;
				}
			}
		}
		for(size_t i=0;i<atom;i++){
			for(size_t j=0;j<3;j++){
				fsout<<now[i][j]<<" ";
			}
			    fsout<<std::endl;
		}
		}
	}
}
