#include "atom.h"
#include <stdio.h>
#include "image.h"
#include <new>
#include <list>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <vector>
double distance(double* a,double* b){
	double s=0.0;
	for(size_t i=0;i<3;i++){
		s=s+(a[i]-b[i])*(a[i]-b[i]);
	}
	return sqrt(s);
}
box::box(){
	virtsize=0;
	p=NULL;
	size=0;
	type=0;
	virtatom=NULL;
	r0=NULL;
	v0=NULL;
	cij=NULL;
	sij=NULL;
	svvij=NULL;
	bvrcut=NULL;
	bvenergy=0.0;
	bvvenergy=0.0;
	ljenergy=0.0;
	epsilonenergy=0.0;
	bvvrcut=NULL;
	vv0=NULL;
	ljrcut=0.0;
	bij=NULL;
	epsilon=NULL;
	stress=NULL;
	stressdft=NULL;
	dftenergy=0.0;
	mdenergy=0.0;
	mdrelevantenergy=0.0;
	dftrelevantenergy=0.0;
	weight=0.0;
}
box::box(atom* inputallatom,
		int t,
		int s,
		double* period,
		double** pairbv_input,
		double** pairbvv_input,
        double** pairlj_input,
        double ljcut
        ){
	/*this paribv_input should be similar to lammps input*/
	/*this paribvv_input should be similar to lammps input*/
	p=new double [3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	std::copy(inputallatom,inputallatom+s,allatom);
	type=t;//specify how many type are in the simulation
	size=s;//specify how many atoms are in the simulation
	/*should have C(n,2) pair*/
	r0=new double* [t];/*equilibrium bond length*/
	v0=new double* [t];/*equilibrium valence zero only diagonal elements are considered*/
	cij=new double* [t];/*power law of bond valence*/
	sij=new double* [t];/*only contains diagonal elements*/
	svvij=new double* [t];/*only contains diagnoal elements*/
	bvrcut=new double* [t];/*cut-off for bond valence*/
	bvvrcut=new double* [t];/*cut-off for bond valence vector*/
	vv0=new double* [t];/*equlibrium bvv0*/
  epsilon = new double* [t];
  bij = new double* [t];
	for(size_t i=0;i<t;i++){
		r0[i]=new double[t];
		v0[i]=new double[t];
		cij[i]=new double[t];
		sij[i]=new double[t];
		svvij[i]=new double[t];
		bvrcut[i]=new double[t];
		bvvrcut[i]=new double[t];
		vv0[i]=new double[t];
    epsilon[i] = new double[t];
    bij[i] = new double[t];
	}
	size_t temp=0;
	double maxcutoff=0.0;
	for(size_t i=0;i<t;i++)
		for(size_t j=i;j<t;j++){
			r0[i][j]=pairbv_input[temp][0];
			r0[j][i]=r0[i][j];
			cij[i][j]=pairbv_input[temp][1];
			cij[j][i]=cij[i][j];
			sij[i][j]=pairbv_input[temp][2];
			sij[j][i]=sij[i][j];
			svvij[i][j]=pairbvv_input[temp][2];
			svvij[j][i]=svvij[i][j];
			v0[i][j]=pairbv_input[temp][3];
			v0[j][i]=v0[i][j];
			vv0[i][j]=pairbvv_input[temp][3];
			vv0[j][i]=vv0[i][j];
			bvrcut[i][j]=pairbv_input[temp][4];
			bvrcut[j][i]=bvrcut[i][j];
			bvvrcut[i][j]=pairbvv_input[temp][4];
			bvvrcut[j][i]=bvvrcut[i][j];
      epsilon[i][j] = pairlj_input[temp][0];
      epsilon[j][i] = epsilon[i][j];
      bij[i][j] = pairlj_input[temp][1];
      bij[j][i] = bij[i][j];
      maxcutoff=maxcutoff > bvrcut[i][j] ? maxcutoff : bvrcut[i][j];
			maxcutoff=maxcutoff > ljrcut ? maxcutoff : ljrcut;
      temp++;
		}
	int virt_size;
	virtatom=imageall(allatom,size,period,maxcutoff,virt_size);
	virtsize=virt_size;
  ljrcut = ljcut;
	stress=new double* [3];
	for(size_t i=0;i<3;i++){
		stress[i]=new double [3];
		for(size_t j=0;j<3;j++){
			stress[i][j]=0.0;
		}
	}
}
/** this input is the same format as control.opt file bond valence optimized parameter looks like the following.
/*
#          r          N          S      alpha          V          R          A          B          C          n       Sbvv        bvv
 1  1    0.00000    5.00000    0.34703    2.00000    2.00000    8.00000    1.00000    2.93318    1.00000   12.00000    6.21787    0.00000 
 1  2    0.00000    5.00000    0.00000    0.00000    0.00000    8.00000    1.00000    2.40135    1.00000   12.00000    0.00000    0.00000 
 1  3    0.00000    5.00000    0.00000    0.00000    0.00000    8.00000    1.00000    1.37294    1.00000   12.00000    0.00000    0.00000 
 1  4    2.35503    8.66909    0.00000    0.00000    0.00000    8.00000    1.00000    1.82651    1.00000   12.00000    0.00000    0.00000 
 2  2    0.00000    5.00000    1.91961    2.00000    2.00000    8.00000    1.00000    2.69873    1.00000   12.00000    5.15915    0.00553 
 2  3    0.00000    5.00000    0.00000    0.00000    0.00000    8.00000    1.00000    2.44973    1.00000   12.00000    0.00000    0.00000 
 2  4    1.91771    7.17080    0.00000    0.00000    0.00000    8.00000    1.00000    1.89220    1.00000   12.00000    0.00000    0.00000 
 3  3    0.00000    5.00000    0.60653    2.00000    4.00000    8.00000    1.00000    3.11268    1.00000   12.00000    0.15963    0.90215 
 3  4    1.86405    4.44881    0.00000    0.00000    0.00000    8.00000    1.00000    1.56763    1.00000   12.00000    0.00000    0.00000 
 4  4    0.00000    5.00000    0.01705    2.00000    2.00000    8.00000    1.00000    1.86348    1.00000   12.00000    5.46101    0.00001 
 * 
 * call function box::updatebvparameter(control::bvvmatrix);
 * */
void box::updatebvparameter(double** input){
	size_t temp=0;	
	for(size_t i=0;i<type;i++)
		for(size_t j=i;j<type;j++){
			r0[i][j]=input[temp][0];
			r0[j][i]=r0[i][j];
			cij[i][j]=input[temp][1];
			cij[j][i]=cij[i][j];
			sij[i][j]=input[temp][2];
			sij[j][i]=sij[i][j];
			svvij[i][j]=input[temp][10];
			svvij[j][i]=svvij[i][j];
			v0[i][j]=input[temp][4];
			v0[j][i]=v0[i][j];
			vv0[i][j]=input[temp][11];
			vv0[j][i]=vv0[i][j];
			bvrcut[i][j]=input[temp][5];
			bvrcut[j][i]=bvrcut[i][j];
			bvvrcut[i][j]=input[temp][5];
			bvvrcut[j][i]=bvvrcut[i][j];
      epsilon[i][j] = epsilon_lj;
      epsilon[j][i] = epsilon[i][j];
      bij[i][j] = input[temp][7];
      bij[j][i] = bij[i][j];
	  	temp++;
		}
}
void box::init(atom* inputallatom,int s,int t,double maxcutoff,double* period,double** input,double dft_energy,double** stress_dft,double w){
	p=new double[3];
	for(size_t i=0;i<3;i++){
		p[i]=period[i];
	}
	allatom=new atom[s];
	std::copy(inputallatom,inputallatom+s,allatom);
	size=s;
	dftenergy=dft_energy;
	weight=w;
	stressdft=new double* [3];
	for(size_t i=0;i<3;i++){
		stressdft[i]=new double [3];
		for(size_t j=0;j<3;j++){
			stressdft[i][j]=stress_dft[i][j];
		}
	}
  type=t;//specify how many type are in the simulation
	size=s;//specify how many atoms are in the simulation
	/*should have C(n,2) pair*/
	r0=new double* [t];/*equilibrium bond length*/
	v0=new double* [t];/*equilibrium valence zero only diagonal elements are considered*/
	cij=new double* [t];/*power law of bond valence*/
	sij=new double* [t];/*only contains diagonal elements*/
	svvij=new double* [t];/*only contains diagnoal elements*/
	bvrcut=new double* [t];/*cut-off for bond valence*/
	bvvrcut=new double* [t];/*cut-off for bond valence vector*/
	vv0=new double* [t];/*equlibrium bvv0*/
    epsilon = new double* [t];
    bij = new double* [t];
	for(size_t i=0;i<t;i++){
		r0[i]=new double[t];
		v0[i]=new double[t];
		cij[i]=new double[t];
		sij[i]=new double[t];
		svvij[i]=new double[t];
		bvrcut[i]=new double[t];
		bvvrcut[i]=new double[t];
		vv0[i]=new double[t];
    epsilon[i] = new double[t];
    bij[i] = new double[t];
	}
	int virt_size;
	virtatom=imageall(allatom,size,period,maxcutoff,virt_size);
	virtsize=virt_size;
  ljrcut = maxcutoff;
	this->updatebvparameter(input);
	this->updatelistbv();
	this->updatelistbvv();
	this->updatelistlj();
}
void box::freezeforce(){
	for(size_t i=0;i<size;i++){
		allatom[i].force[0]=0.0;
		allatom[i].force[1]=0.0;
		allatom[i].force[2]=0.0;
	}
}
void box::computeAll(){//zhenbang
    freezeforce();
    computebv();
    computebvv();
    //computestress();
    computelj();
    computelong(1e-10);
		mdenergy=bvenergy+bvvenergy+ljenergy+epsilonenergy;
}

void box::printnei(int i){
   for(std::list<int>::iterator a=allatom[i].neibv.begin();a!=allatom[i].neibv.end();a++){
      std::cout<<" "<<*a;
   }
}
/*finished computing bond valence force*/
/*end define the bond-valence energy*/
void box::printlj(){
	for(size_t i=0;i<size;i++){
		std::cout<<allatom[i].force[0]<<" "<<allatom[i].force[1]<<" "<<allatom[i].force[2]<<std::endl;
	}
}
void box::printinfo(){
	for(size_t i=0;i<size;i++){
		for(size_t j=0;j<3;j++){
			std::cout<<allatom[i].position[j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<p[0]<<" "<<p[1]<<" "<<p[2]<<std::endl;
	for(size_t i=0;i<size;i++){
		for(size_t j=0;j<3;j++){
			std::cout<<allatom[i].dftforce[j]<<" ";
		}
		std::cout<<std::endl;
	}
	std::cout<<dftenergy<<std::endl;
	for(size_t i=0;i<3;i++){
		for(size_t j=0;j<3;j++){
			std::cout<<stressdft[i][j]<<" ";
		}
	std::cout<<std::endl;
	}
}
void box::settype(std::vector<int>& t){
	size_t temp=0;
	for(size_t i=0;i<t.size();i++)
		for(size_t j=0;j<t[i];j++){
			allatom[i].type=i;
			temp++;
		}
}
