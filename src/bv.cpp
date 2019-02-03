#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
#include <math.h>
void box::updatelistbv(){
	double temp;
	double paircut;
	for(size_t i=0;i<size;i++){
		allatom[i].neibv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
            paircut=bvrcut[allatom[i].type][virtatom[j].type];
            if(temp<paircut && temp>0.0000001){
				allatom[i].neibv.push_back(j);
			}
		}
	}
}
/*starting a light version of bond valence with computing the force and bond valence energy*/
void box::computebv(){
	/*bond valence parameters.
	 *typeone typetwo r0 Nij S V0 rcut
	 * 1     1        #  #  # #  #
	 * 1     1        #  #  # #  #
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 * 1     1        ###############
	 */
	double delx,dely,delz,rsq,recip,r,s;
	bvenergy=0.00;
	double* fp=new double[size];
	for(size_t i=0;i<size;i++){
		allatom[i].s0=0;
        //        std::cout<<allatom[i].neibv.size()<<std::endl;
		for(std::list<int>::iterator j=allatom[i].neibv.begin();j!=allatom[i].neibv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			allatom[i].s0+=pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type]);
		}
		s=allatom[i].s0-v0[allatom[i].type][allatom[i].type];
    bvenergy=sij[allatom[i].type][allatom[i].type]*(s*s)+bvenergy;
		fp[i]=2*sij[allatom[i].type][allatom[i].type]*s;
	}
	/*finished computing energy and started to compute force*/
        double Aij=0.0;
	for(size_t i=0;i<size;i++){
		for(std::list<int>::iterator j=allatom[i].neibv.begin();j!=allatom[i].neibv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			Aij=cij[allatom[i].type][virtatom[*j].type]*pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type])/r;
			allatom[i].force[0]+=(fp[i]+fp[*j%size])*Aij*delx/r;
			allatom[i].force[1]+=(fp[i]+fp[*j%size])*Aij*dely/r;
			allatom[i].force[2]+=(fp[i]+fp[*j%size])*Aij*delz/r;
		}
    }
}
void box::updatebv(double** pairbv_input){
			size_t temp=0;
				for(size_t i=0;i<type;i++)
							for(size_t j=i;j<type;j++){
									r0[i][j]=pairbv_input[temp][0];
									r0[j][i]=r0[i][j];
								  cij[i][j]=pairbv_input[temp][1];
								  r0[j][i]=r0[i][j];
								  sij[i][j]=pairbv_input[temp][2];
								  sij[j][i]=svvij[i][j];
								  v0[i][j]=pairbv_input[temp][3];
								  v0[j][i]=v0[i][j];
								 bvrcut[i][j]=pairbv_input[temp][4];
						     bvrcut[j][i]=bvvrcut[i][j];
			           temp++;
							}
}
