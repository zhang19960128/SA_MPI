#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
#include <math.h>
void box::updatelistbvv(){
	double temp;
	double paircut;
	for(size_t i=0;i<size;i++){
		allatom[i].neibvv.clear();
		for(size_t j=0;j<virtsize;j++){
			temp=distance(allatom[i].position,virtatom[j].position);
            paircut=bvvrcut[allatom[i].type][virtatom[j].type];
            if(temp<paircut && temp>0.0000001){
				allatom[i].neibvv.push_back(j);
			}
		}
	}
}
void box::updatebvv(double** pairbvv_input){
	size_t temp=0;
	for(size_t i=0;i<type;i++)
		for(size_t j=i;j<type;j++){
			r0[i][j]=pairbvv_input[temp][0];
                        r0[j][i]=r0[i][j];
			cij[i][j]=pairbvv_input[temp][1];
                        r0[j][i]=r0[i][j];
			svvij[i][j]=pairbvv_input[temp][2];
                        svvij[j][i]=svvij[i][j];
			vv0[i][j]=pairbvv_input[temp][3];
                        vv0[j][i]=v0[i][j];
			bvvrcut[i][j]=pairbvv_input[temp][4];
                        bvvrcut[j][i]=bvvrcut[i][j];
			temp++;
		}
}
void box::computebvv(){
   /* bond valence vector parameters
    * typeone typetwo r0 Nij S V00 rcut
    *
    *
    *
    *
    *
    *
    */
  double delx,dely,delz,rsq,r,s,ss;
	bvvenergy=0.00;
	double* fp=new double[size];
        double temp;
				double recip;
				double recip2;
				double Aij;
				double Eij;
				double fx;
				double fy;
				double fz;
	double* Dix=new double [size];//store some temperory variables
	double* Diy=new double [size];//
	double* Diz=new double [size];//
	int itype,jtype;
	for(size_t i=0;i<size;i++){
		allatom[i].s0x=0.0;
		allatom[i].s0y=0.0;
		allatom[i].s0z=0.0;
		for(std::list<int>::iterator j=allatom[i].neibvv.begin();j!=allatom[i].neibvv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			r=sqrt(rsq);
			recip=1.0/r;
			temp=pow(r0[allatom[i].type][virtatom[*j].type]/r,cij[allatom[i].type][virtatom[*j].type]);
			allatom[i].s0x+=temp*delx/r;
		//	if(i==0){
		//		std::cout<<i<<" "<<std::setprecision(9)<<allatom[i].s0x<<std::endl;;
		//	}
			allatom[i].s0y+=temp*dely/r;
			allatom[i].s0z+=temp*delz/r;
		}
		s=allatom[i].s0x*allatom[i].s0x+allatom[i].s0y*allatom[i].s0y+allatom[i].s0z*allatom[i].s0z-vv0[allatom[i].type][allatom[i].type]*vv0[allatom[i].type][allatom[i].type];
		ss=s*s;
		bvvenergy=bvvenergy+svvij[allatom[i].type][allatom[i].type]*ss;
		Dix[i]=svvij[allatom[i].type][allatom[i].type]*2*2*allatom[i].s0x*s;
		Diy[i]=svvij[allatom[i].type][allatom[i].type]*2*2*allatom[i].s0y*s;
		Diz[i]=svvij[allatom[i].type][allatom[i].type]*2*2*allatom[i].s0z*s;
	//	std::cout<<std::fixed<<std::setprecision(10)<<"\t"<<std::setw(15)<<s<<" "<<std::setw(15)<<Dix[i]<<" "<<std::setw(15)<<Diy[i]<<" "<<std::setw(15)<<Diz[i]<<std::endl;
	}
//	std::cout<<"the total energy is "<<std::setprecision(10)<<bvenergy+bvvenergy<<std::endl;
	/*finished calculating the bond valence vector energy*/
  /*started to compute the */
	for(size_t i=0;i<size;i++){
		for(std::list<int>::iterator j=allatom[i].neibvv.begin();j!=allatom[i].neibvv.end();j++){
			delx=allatom[i].position[0]-virtatom[*j].position[0];
			dely=allatom[i].position[1]-virtatom[*j].position[1];
			delz=allatom[i].position[2]-virtatom[*j].position[2];
			rsq=delx*delx+dely*dely+delz*delz;
			itype=allatom[i].type;
			jtype=virtatom[*j].type;
			r=sqrt(rsq);
			recip=1.0/r;
			recip2=recip*recip;
			Aij=pow(r0[itype][jtype]*recip,cij[itype][jtype])*recip;
			Eij=(cij[itype][jtype]+1.0)*recip2;
        fx = (Dix[*j%size]-Dix[i])*Aij
             + (Dix[i]-Dix[*j%size])*Eij*delx*delx*Aij
             + (Diy[i]-Diy[*j%size])*Eij*delx*dely*Aij
             + (Diz[i]-Diz[*j%size])*Eij*delx*delz*Aij;
        fy = (Diy[*j%size]-Diy[i])*Aij
             + (Diy[i]-Diy[*j%size])*Eij*dely*dely*Aij
             + (Diz[i]-Diz[*j%size])*Eij*dely*delz*Aij
             + (Dix[i]-Dix[*j%size])*Eij*dely*delx*Aij;
        fz = (Diz[*j%size]-Diz[i])*Aij
             + (Diz[i]-Diz[*j%size])*Eij*delz*delz*Aij
             + (Dix[i]-Dix[*j%size])*Eij*delz*delx*Aij
             + (Diy[i]-Diy[*j%size])*Eij*delz*dely*Aij;
				allatom[i].force[0]+=fx;
				allatom[i].force[1]+=fy;
				allatom[i].force[2]+=fz;
		}
	//	std::cout<<std::fixed<<std::setprecision(15)<<std::setw(15)<<allatom[i].force[0]<<" "<<std::setw(15)<<allatom[i].force[1]<<" "<<std::setw(15)<<allatom[i].force[2]<<std::endl;
	}
	delete [] Dix;
	delete [] Diy;
	delete [] Diz;
}
