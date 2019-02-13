#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <iomanip>
#include "readpara.h"
#define pi 3.14159265359
/*Perform the Ewald summation. The lattice paramters should be given*/
#define EWALD_F   1.12837917//actually this is 2/sqrt(3)
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429
/*Lattice should be a flattened array*/
double rms(double g_ewald,int km,double prd,int natoms,double q2){
	double value=2.0*q2*g_ewald/prd*sqrt(1.0/pi/km/natoms)*exp(-1*pi*pi*km*km/(g_ewald*g_ewald*prd*prd));
	return value;
}
void box::computelong(double accuracy_relative){
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double epsil = 0.0055263885;
		double coul_prefactor=180.9512801302711;
		double ewald_alpha;
		/*e/(epsilon0*A)=180.951  ev when use U=k*q1*q2/r*/
    double volume=p[0]*p[1]*p[2];
    double delx,dely,delz,rsq,r;
    double root2pi=sqrt(2*pi);
		double root2=sqrt(2);
		double rootpi=sqrt(pi);
    double* xall=new double[size];
    double* yall=new double [size];
    double* zall=new double [size];
		double* fx=new double [size];
		double* fy=new double [size];
		double* fz=new double [size];
		/**
		 * pass the charge from control::charge to here
		 */
		/******************************************************/
		double* chargetype=new double [species::spe.size()];
		for(size_t i=0;i<species::spe.size();i++){
			chargetype[i]=control::charge[i];
		}
		/******************************************************/
		double chargei,chargej,temp,temp2,r3,erfc_interpolate,expm2,grij,t;//erfc_exact; use to debug when compare different erfc function.
  /*determin the value of g_ewald according to the accuracy you want*/
	double q2=0.0;
	double accuracy=accuracy_relative*coul_prefactor/4.0/pi;
	for(size_t i=0;i<size;i++){
		q2=chargetype[allatom[i].type]*chargetype[allatom[i].type]+q2;
	}
	q2=q2*coul_prefactor/4/pi;
	ewald_alpha=accuracy*sqrt(size*ljrcut*p[0]*p[1]*p[2])/2.0/q2;
	if(ewald_alpha >=1.0){
	   ewald_alpha=(1.35-0.15*log(accuracy))/ljrcut;
	}
	else{
		ewald_alpha=sqrt(-log(ewald_alpha))/ljrcut;
	}
	int gmax=1;
	double err=rms(ewald_alpha,gmax,p[0],size,q2);
    while (err > accuracy) {
      gmax++;
      err = rms(ewald_alpha,gmax,p[0],size,q2);
    }
    while (err > accuracy) {
      gmax++;
      err = rms(ewald_alpha,gmax,p[1],size,q2);
    }  
		while (err > accuracy) {
      gmax++;
      err = rms(ewald_alpha,gmax,p[1],size,q2);
    }
		ewald_alpha=0.684653;
		gmax=9;
		double sigma=1.0/sqrt(2)/ewald_alpha;
		/****finish esitimate the g_ewald and kmax******/
		for(size_t i=0;i<size;i++){
       xall[i]=allatom[i].position[0];
       yall[i]=allatom[i].position[1];
       zall[i]=allatom[i].position[2];
			 fx[i]=0.00;
			 fy[i]=0.00;
			 fz[i]=0.00;
       chargei=chargetype[allatom[i].type];
			 //std::cout<<"I got "<<allatom[i].neilj.size()<<"neighbors"<<std::endl;
       for(std::list<int>::iterator j=allatom[i].neilj.begin();j!=allatom[i].neilj.end();j++){
         chargej=chargetype[virtatom[*j].type];
				 delx=allatom[i].position[0]-virtatom[*j].position[0];
         dely=allatom[i].position[1]-virtatom[*j].position[1];
         delz=allatom[i].position[2]-virtatom[*j].position[2];
         rsq=delx*delx+dely*dely+delz*delz;
         r=sqrt(rsq);
				 r3=r*rsq;
				// temp=erfc(r/root2/sigma);
        // ShortRange+=1/epsil/4/pi*chargei*chargej/r*temp;
				 //temp2=1.0/8.0/epsil/pi*chargei*chargej/r3*erfc(r/root2/sigma); this is the standard erfc function but now I want to use the interpolation method to accelerate.
				 grij=ewald_alpha*r;
				 expm2=exp(-1.0*grij*grij);//=exp(-1*rsq/2/sigma^2)
				 t= 1.0 / (1.0 + EWALD_P*grij);
				 erfc_interpolate = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
				 //std::cout<<"the difference is: "<<erfc_interpolate-erfc(r/root2/sigma)<<std::endl;
				 ShortRange+=1.0/epsil/4.0/pi*chargei*chargej/r*erfc_interpolate;
				 temp=1.0/4.0/pi/epsil*chargei*chargej/r3*(EWALD_F*expm2*grij+erfc_interpolate);
				 fx[i]=fx[i]+temp*delx;
				 fy[i]=fy[i]+temp*dely;
				 fz[i]=fz[i]+temp*delz;
			 }
			 //std::cout<<i<<" "<<fx[i]<<" "<<fy[i]<<" "<<fz[i]<<std::endl;
    }
    ShortRange=ShortRange/2.0;
    for(size_t i=0;i<size;i++){
      selfe=selfe-1/root2pi/sigma/4/pi/epsil*(chargetype[allatom[i].type])*chargetype[allatom[i].type];
    }
		//std::cout<<"the self energy is: "<<selfe<<std::endl;
    /*calculate the longrange energy*/
    /*exp(i*k_v*r)
     *=exp(i*k_x*r_x)*exp(i*k_y*r_y)*exp(i*k_z*r_z)
     *=cos(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*sin(k_y*r_y)*cos(k_z*r_z)-sin(k_x*r_x)*cos(k_y*r_y)*cos(k_z*r_z)-cos(k_x*r_x)*sin(k_y*r_y)*sin(k_z*r_z)
     *+I(********************************************************) similar term just expand it
     *
     * */
    /*define array for cos(k_x*r_x),cos(k_y*r_y),cos(k_z*r_z),sin(k_x*r_x),sin(k_y*r_y),sin(k_z*r_z)*/
    double** cskxrx=new double* [gmax+1];
    double** cskyry=new double* [gmax+1];
    double** cskzrz=new double* [gmax+1];
    double** snkxrx=new double* [gmax+1];
    double** snkyry=new double* [gmax+1];
    double** snkzrz=new double* [gmax+1];
    for(size_t i=0;i<=gmax;i++){
      cskxrx[i]=new double [size];
      cskyry[i]=new double [size];
      cskzrz[i]=new double [size];
      snkxrx[i]=new double [size];
      snkyry[i]=new double [size];
      snkzrz[i]=new double [size];
    }
    for(int i=0;i<=gmax;i++)
       for(size_t j=0;j<size;j++){
          cskxrx[i][j]=cos(2*pi/p[0]*i*xall[j]);
          cskyry[i][j]=cos(2*pi/p[1]*i*yall[j]);
          cskzrz[i][j]=cos(2*pi/p[2]*i*zall[j]);
          snkxrx[i][j]=sin(2*pi/p[0]*i*xall[j]);
          snkyry[i][j]=sin(2*pi/p[1]*i*yall[j]);
          snkzrz[i][j]=sin(2*pi/p[2]*i*zall[j]);
			//		std::cout<<cskxrx[i][j]<<" "<<cskyry[i][j]<<" "<<cskzrz[i][j]<<std::endl;
				  //std::cout<<xall[j]<<" "<<yall[j]<<" "<<zall[j]<<std::endl;
       }
    double skre,skim,skmodsq,ksq,snkr,cskr,kr;
    for(int h=-1*gmax;h<=gmax;h++)
       for(int k=-1*gmax;k<=gmax;k++)
          for(int l=-1*gmax;l<=gmax;l++){
             skre=0.0;
             skim=0.0;
	     if(h==0&&k==0&&l==0){
		  continue;
	       }
             for(size_t i=0;i<size;i++){
	        /*finish computing real part*/
							 chargei=chargetype[allatom[i].type];
                skre=skre+(cskxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i])*chargei;
	        			skre=skre-chargei*snkxrx[abs(h)][i]*snkyry[abs(k)][i]*cskzrz[abs(l)][i]*(h>0 ? 1:-1)*(k>0 ? 1:-1);
								skre=skre-chargei*snkxrx[abs(h)][i]*cskyry[abs(k)][i]*snkzrz[abs(l)][i]*(h>0 ? 1:-1)*(l>0 ? 1:-1);
                skre=skre-chargei*cskxrx[abs(h)][i]*snkyry[abs(k)][i]*snkzrz[abs(l)][i]*(k>0 ? 1:-1)*(l>0 ? 1:-1);
								/*start to compute the imaginary part*/
	        			skim=skim+chargei*cskxrx[abs(h)][i]*cskyry[abs(k)][i]*snkzrz[abs(l)][i]*(l>0 ? 1:-1);
	        			skim=skim+chargei*snkxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i]*(h>0 ? 1:-1);
	        			skim=skim+chargei*cskxrx[abs(h)][i]*snkyry[abs(k)][i]*cskzrz[abs(l)][i]*(k>0 ? 1:-1);
	        			skim=skim-chargei*snkxrx[abs(h)][i]*snkyry[abs(k)][i]*snkzrz[abs(l)][i]*(h>0 ? 1:-1)*(k>0 ? 1:-1)*(l>0 ? 1:-1);
						 }
          ksq=h*2*pi/p[0]*h*2*pi/p[0]+k*2*pi/p[1]*k*2*pi/p[1]+l*2*pi/p[2]*l*2*pi/p[2];
		      skmodsq=skre*skre+skim*skim;
		      temp=exp(-1*sigma*sigma*ksq/2)/ksq;
		      LongRange=LongRange+1/volume/2/epsil*temp*skmodsq;
						 for(size_t i=0;i<size;i++){
		//std::cout<<"the real part is: "<<skre<<" the imaginary part is: "<<skim<<std::endl;
		 			chargei=chargetype[allatom[i].type];
		 			kr=h*2*pi/p[0]*xall[i]+k*2*pi/p[1]*yall[i]+l*2*pi/p[2]*zall[i];
					snkr=sin(kr);
					//snkr=(snkxrx[abs(h)][i]*cskyry[abs(k)][i])*cskzrz[abs(l)][i]*(h>0?1:-1);
					//snkr=(snkxrx[abs(h)][i]*cskyry[abs(k)][i]*(h>0?1:-1)+cskxrx[abs(h)][i]*snkyry[abs(k)][i]*(k>0?1:-1))*cskzrz[abs(l)][i];
					//snkr=snkr
					/*
					 * the full formula require this, but due to symmetry, the net sum of snkzrz is basically zeros, so this is useless
					snkr=snkr+(cskxrx[abs(h)][i]*cskyry[abs(k)][i]-snkxrx[abs(h)][i]*snkyry[abs(k)][i]*(h>0?1:-1)*(k>0?1:-1))*snkzrz[abs(l)][i];
		 			*/
					cskr=cos(kr);
				  //	cskr=cskxrx[abs(h)][i]*cskyry[abs(k)][i]*cskzrz[abs(l)][i];
					//  cskr=(cskxrx[abs(h)][i]*cskyry[abs(k)][i]-snkxrx[abs(h)][i]*snkyry[abs(k)][i]*(h>0?1:-1)*(k>0?1:-1))*cskzrz[abs(l)][i];
					//	cskr=cskr-(snkxrx[abs(h)][i]*cskyry[abs(k)][i]*(h>0?1:-1)+cskxrx[abs(h)][i]*snkyry[abs(k)][i]*(k>0?1:-1))*snkzrz[abs(l)][i];
					fx[i]=fx[i]+1/volume/2/epsil*temp*chargei*2*pi/p[0]*h*2*(snkr*skre-cskr*skim);
					    fy[i]=fy[i]+1/volume/2/epsil*temp*chargei*2*pi/p[1]*k*2*(snkr*skre-cskr*skim);
					    fz[i]=fz[i]+1/volume/2/epsil*temp*chargei*2*pi/p[2]*l*2*(snkr*skre-cskr*skim);
						 }
						 }
		epsilonenergy=selfe+ShortRange+LongRange;
		for(size_t i=0;i<size;i++){
			allatom[i].force[0]+=fx[i];
			allatom[i].force[1]+=fy[i];
			allatom[i].force[2]+=fz[i];
		}
		/*this is the final step*/
    for(size_t i=0;i<gmax;i++){
      delete [] cskxrx[i];
      delete [] cskyry[i];
      delete [] cskzrz[i];
      delete [] snkxrx[i];
      delete [] snkyry[i];
      delete [] snkzrz[i];
    }
    delete [] cskxrx;
    delete [] cskyry;
    delete [] cskzrz;
    delete [] snkxrx;
    delete [] snkyry;
    delete [] snkzrz;
    delete [] xall;
    delete [] yall;
    delete [] zall;
		delete [] fx;
		delete [] fy;
		delete [] fz;
		delete [] chargetype;
    /*end memory allocation*/
}
