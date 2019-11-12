#include "atom.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <list>
#include <iomanip>
#include "readpara.h"
#define PI 3.14159265359
/*Perform the Ewald summation. The lattice paramters should be given*/
#define EWALD_F   1.12837917//actually this is 2/sqrt(3)
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429
#define _TABULATE_SINE
#define _TABULATE_ERFC
/* real space positions and charges */
double alpha,rcut,kcut; 
/* H matrix */
double h11,h12,h13,h21,h22,h23,h31,h32,h33; 
/* G matrix: H*G^T = 2\PI I */
double g11,g12,g13,g21,g22,g23,g31,g32,g33,volume;
/* reciprocal lattices to be be summed */
int max_g1,max_g2,max_g3;
/* is the (x,y,z) component of the first cell edge, etc. */ 
void init_cell (double** H)
{ 
    /* Because H is assumed to be coming from Fortran */
    /* we have to do a little conversion */
    h11 = H[0][0];  /* 1st element */
    h12 = H[0][1];  /* 4th element */
    h13 = H[0][2];  /* 7th element */
    h21 = H[1][0];  /* 2nd element */
    h22 = H[1][1];  /* 5th element */
    h23 = H[1][2];  /* 8th element */
    h31 = H[2][0];  /* 3rd element */
    h32 = H[2][1];  /* 6th element */
    h33 = H[2][2];  /* 9th element */
    /* Get the reciprocal lattice vectors */
    /*calculated using Mathematica*/
    g11=(-(h23*h32) + h22*h33)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 -h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g12=(h23*h31 - h21*h33)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g13=(-(h22*h31) + h21*h32)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g21=(h13*h32 - h12*h33)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g22=(-(h13*h31) + h11*h33)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g23=(h12*h31 - h11*h32)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 -h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g31=(-(h13*h22) + h12*h23)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g32=(h13*h21 - h11*h23)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    g33=(-(h12*h21) + h11*h22)/(-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33);
    volume=-(h13*h22*h31) + h12*h23*h31 + h13*h21*h32 - h11*h23*h32 - h12*h21*h33 + h11*h22*h33;
    return;
}  /* end init_cell() */
double compute_ksq(int i,int j,int k){
  return (i*g11+j*g21+k*g31)*(i*g11+j*g21+k*g31)+(i*g12+j*g22+k*g32)*(i*g12+j*g22+k*g32)+(i*g13+j*g23+k*g33)*(i*g13+j*g23+k*g33);
}
int fetch(int i,int g){
  return i-g;
}
/*Lattice should be a flattened array*/
double rms(double g_ewald,int km,double prd,int size,double q2){
	double value=2.0*q2*g_ewald/prd*sqrt(1.0/PI/km/size)*exp(-1*PI*PI*km*km/(g_ewald*g_ewald*prd*prd));
	return value;
}
void box::computelong(double accuracy_relative){
    double duration;
    double ShortRange = 0;
    double LongRange = 0;
    double selfe = 0;
    double epsil = 0.0055263885;
		double coul_prefactor=180.9512801302711;
		double ewald_alpha;
		/*e/(epsilon0*A)=180.951  ev when use U=k*q1*q2/r*/
    double delx,dely,delz,rsq,r;
    double root2PI=sqrt(2*PI);
		double root2=sqrt(2);
		double rootPI=sqrt(PI);
    double* xall=new double[size];
    double* yall=new double [size];
    double* zall=new double [size];
    double* s1=new double [size];
    double* s2=new double [size];
    double* s3=new double [size];
		double* fx=new double [size];
		double* fy=new double [size];
		double* fz=new double [size];
    double* q=new double [size];
		/**
		 * pass the charge from control::charge to here
		 */
		/******************************************************/
		double* chargetype=new double [species::spe.size()];
		for(int  i=0;i<species::spe.size();i++){
			chargetype[i]=control::charge[i];
		}
		/******************************************************/
		double chargei,chargej,temp,temp2,r3,erfc_interpolate,expm2,grij,t;//erfc_exact; use to debug when compare different erfc function.
    init_cell(p);
  /*determin the value of g_ewald according to the accuracy you want*/
	double q2=0.0;
	double accuracy=accuracy_relative*coul_prefactor/4.0/PI;
	for(int  i=0;i<size;i++){
		q2=chargetype[allatom[i].type]*chargetype[allatom[i].type]+q2;
	}
	q2=q2*coul_prefactor/4/PI;
	ewald_alpha=accuracy*sqrt(size*ljrcut*p[0][0]*p[1][1]*p[2][2])/2.0/q2;
	if(ewald_alpha >=1.0){
	   ewald_alpha=(1.35-0.15*log(accuracy))/ljrcut;
	}
	else{
		ewald_alpha=sqrt(-log(ewald_alpha))/ljrcut;
	}
 int   max_g1 = 1;
 int   max_g2 = 1;
 int   max_g3 = 1;
 double err;
    err = rms(ewald_alpha,max_g1,p[0][0],size,q2);
    while (err > accuracy) {
      max_g1++;
      err = rms(ewald_alpha,max_g1,p[0][0],size,q2);
    }

    err = rms(ewald_alpha,max_g2,p[1][1],size,q2);
    while (err > accuracy) {
      max_g2++;
      err = rms(ewald_alpha,max_g2,p[1][1],size,q2);
    }

    err = rms(ewald_alpha,max_g3,p[2][2],size,q2);
    while (err > accuracy) {
      max_g3++;
      err = rms(ewald_alpha,max_g3,p[2][2],size,q2);
    }
  ewald_alpha=0.684653;
	double sigma=1.0/sqrt(2)/ewald_alpha;
  double rij;
		/****finish esitimate the g_ewald and kmax******/
		for(int  i=0;i<size;i++){
       xall[i]=allatom[i].position[0];
       yall[i]=allatom[i].position[1];
       zall[i]=allatom[i].position[2];
       s1[i]=allatom[i].crystal_position[0];
       s2[i]=allatom[i].crystal_position[1];
       s3[i]=allatom[i].crystal_position[2];
			 fx[i]=0.00;
			 fy[i]=0.00;
			 fz[i]=0.00;
       chargei=chargetype[allatom[i].type];
       q[i]=chargei;
       for(std::list<int>::iterator j=allatom[i].neilj.begin();j!=allatom[i].neilj.end();j++){
         chargej=chargetype[virtatom[*j].type];
				 delx=allatom[i].position[0]-virtatom[*j].position[0];
         dely=allatom[i].position[1]-virtatom[*j].position[1];
         delz=allatom[i].position[2]-virtatom[*j].position[2];
         rsq=delx*delx+dely*dely+delz*delz;
         r=sqrt(rsq);
				 r3=r*rsq;
				 grij=ewald_alpha*r;
				 expm2=exp(-1.0*grij*grij);//=exp(-1*rsq/2/sigma^2)
				 t= 1.0 / (1.0 + EWALD_P*grij);
				 erfc_interpolate = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
				 ShortRange+=1.0/epsil/4.0/PI*chargei*chargej/r*erfc_interpolate;
				 temp=1.0/4.0/PI/epsil*chargei*chargej/r3*(EWALD_F*expm2*grij+erfc_interpolate);
				 fx[i]=fx[i]+temp*delx;
				 fy[i]=fy[i]+temp*dely;
				 fz[i]=fz[i]+temp*delz;
			 }
    }
    ShortRange=ShortRange/2.0;
    for(int  i=0;i<size;i++){
      selfe=selfe-1/root2PI/sigma/4/PI/epsil*(chargetype[allatom[i].type])*chargetype[allatom[i].type];
    }
    /*sin(k dot r),cos( k dot r )*/
    double**** sinkdotr=new double*** [2*max_g1+1];
    double**** coskdotr=new double*** [2*max_g1+1];
    for(size_t i=0;i<2*max_g1+1;i++){
      sinkdotr[i]=new double** [2*max_g2+1];
      coskdotr[i]=new double** [2*max_g2+1];
      for(size_t j=0;j<2*max_g2+1;j++){
        sinkdotr[i][j]=new double* [2*max_g3+1];
        coskdotr[i][j]=new double* [2*max_g3+1];
        for(size_t k=0;k<2*max_g3+1;k++){
          sinkdotr[i][j][k]=new double [size];
          coskdotr[i][j][k]=new double [size];
        }
      }
    }
    double kdotr=0.0;
    for(int i=-max_g1;i<=max_g1;i++){
      for(int j=-max_g2;j<=max_g2;j++){
        for(int k=-max_g3;k<=max_g3;k++){
          for(int m=0;m<size;m++){
          kdotr=2*PI*(s1[m]*i+s2[m]*j+s3[m]*k);/*crytal coordinates and Fourier space are orthogonal*/
          sinkdotr[i+max_g1][j+max_g2][k+max_g3][m]=sin(kdotr);
          coskdotr[i+max_g1][j+max_g2][k+max_g3][m]=cos(kdotr);
          }
        }
      }
    }
    LongRange=0.0;
    double ksq=0.0;
    double prefactor=1/2/volume/epsil;
       /*refer to the notes on OneNote Simulated Annealing,m is the index m in the node, t is the index j in the node*/
    double sinkdotrm,sinkdotrj,coskdotrm,coskdotrj;
    double SQ_real,SQ_imag;//Define SQ=Sigma qi x exp(ikr);
      /*this loop calculate all the forces on each atom and simutaneousl calculate the long Range energy*/
         for(int i=-max_g1;i<=max_g1;i++){
          for(int j=-max_g2;j<=max_g2;j++){
            for(int k=-max_g3;k<=max_g3;k++){
              if(i==0 && j==0 && k==0) continue;
               SQ_real=0.0;
               SQ_imag=0.0;
              ksq=compute_ksq(i,j,k);
              ksq=(2*PI)*(2*PI)*ksq;
               for(size_t m=0;m<size;m++){
                 SQ_real=SQ_real+q[m]*coskdotr[i+max_g1][j+max_g2][k+max_g3][m];
                 SQ_imag=SQ_imag+q[m]*sinkdotr[i+max_g1][j+max_g2][k+max_g3][m];
               }
              prefactor=1.0/2/volume/epsil*exp(-ksq*sigma*sigma/2.0)/ksq;
              LongRange=LongRange+prefactor*(SQ_real*SQ_real+SQ_imag*SQ_imag);
               for(size_t m=0;m<size;m++){
              sinkdotrm=sinkdotr[i+max_g1][j+max_g2][k+max_g3][m];
              coskdotrm=coskdotr[i+max_g1][j+max_g2][k+max_g3][m];
              fx[m]=fx[m]-prefactor*2*PI*(i*g11+j*g21+k*g31)*2*q[m]*(coskdotrm*SQ_imag-sinkdotrm*SQ_real);
              fy[m]=fy[m]-prefactor*2*PI*(i*g12+j*g22+k*g32)*2*q[m]*(coskdotrm*SQ_imag-sinkdotrm*SQ_real);
              fz[m]=fz[m]-prefactor*2*PI*(i*g13+j*g23+k*g33)*2*q[m]*(coskdotrm*SQ_imag-sinkdotrm*SQ_real);
            }
           }
        }
      }
    /*An alternative way to calculate it for your convinience*/
		epsilonenergy=selfe+ShortRange+LongRange;
/*delete all the pointers*/
   for(size_t i=0;i<size;i++){
    allatom[i].force[0]+=fx[i];
    allatom[i].force[1]+=fy[i];
    allatom[i].force[2]+=fz[i];
   }
     for(size_t i=0;i<2*max_g1+1;i++){
      for(size_t j=0;j<2*max_g2+1;j++){
        for(size_t k=0;k<2*max_g3+1;k++){
          delete [] coskdotr[i][j][k];
          delete [] sinkdotr[i][j][k];
        }
        delete [] coskdotr[i][j];
        delete [] sinkdotr[i][j];
      }
      delete [] coskdotr[i];
      delete [] sinkdotr[i];
    }
    delete [] coskdotr;
    delete [] sinkdotr;
    delete [] xall;
    delete [] yall;
    delete [] zall;
		delete [] fx;
		delete [] fy;
		delete [] fz;
		delete [] chargetype;
    delete [] s1;
    delete [] s2;
    delete [] s3;
    delete [] q;
    /*end memory allocation*/
}
