#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI acos(-1.0)

/* number of particles */
static int NP;
/* real space positions and charges */
static double *x,*y,*z,*q;
/* Ewald summation parameters */
static double my_accuracy,alpha,rcut,kcut; 
/* H matrix */
static double h11,h12,h13,h21,h22,h23,h31,h32,h33; 
/* G matrix: H*G^T = 2\pi I */
static double g11,g12,g13,g21,g22,g23,g31,g32,g33,volume;
/* Distance from origin to three unit    */
/* cell surfaces and reciprocal surfaces */
static double rd1,rd2,rd3,kd1,kd2,kd3;
/* real space lattices to be be summed */
static int max_n1,max_n2,max_n3,max_xvec,num_xvec,*nvec;
static double *xvec;
/* reciprocal lattices to be be summed */
static int max_g1,max_g2,max_g3,max_gvec,num_gvec,*gvec; 
#ifdef _TABULATE_ERFC
static double derr,*err,*cli;
#endif
#ifdef _TABULATE_SINE
static double *sinn,*cosn;
#endif

/* Our convention is that (h11,h12,h13) (in f77 index)   */
/* is the (x,y,z) component of the first cell edge, etc. */ 
void init_cell (double H[3][3])
{ 
    /* Because H is assumed to be coming from Fortran */
    /* we have to do a little conversion */
    h11 = H[0][0];  /* 1st element */
    h12 = H[1][0];  /* 4th element */
    h13 = H[2][0];  /* 7th element */
    h21 = H[0][1];  /* 2nd element */
    h22 = H[1][1];  /* 5th element */
    h23 = H[2][1];  /* 8th element */
    h31 = H[0][2];  /* 3rd element */
    h32 = H[1][2];  /* 6th element */
    h33 = H[2][2];  /* 9th element */
    /* Get the reciprocal lattice vectors */
    g11 = h22*h33 - h23*h32;
    g22 = h33*h11 - h31*h13;
    g33 = h11*h22 - h12*h21;
    g12 = h23*h31 - h21*h33;
    g23 = h31*h12 - h32*h11;
    g31 = h12*h23 - h13*h22;
    g13 = h21*h32 - h31*h22;
    g21 = h32*h13 - h12*h33;
    g32 = h13*h21 - h23*h11;
    volume = h11*g11 + h12*g12 + h13*g13;
    /* the shortest distance to respective surfaces */
    rd1 = fabs(volume)/sqrt(g11*g11+g12*g12+g13*g13);
    rd2 = fabs(volume)/sqrt(g21*g21+g22*g22+g23*g23);
    rd3 = fabs(volume)/sqrt(g31*g31+g32*g32+g33*g33);
    /* reciprocal lattice vectors */
    g11 *= 2*PI/volume;
    g12 *= 2*PI/volume;
    g13 *= 2*PI/volume;
    g21 *= 2*PI/volume;
    g22 *= 2*PI/volume;
    g23 *= 2*PI/volume;
    g31 *= 2*PI/volume;
    g32 *= 2*PI/volume;
    g33 *= 2*PI/volume;
    /* the shortest distance to respective reciprocal surfaces */
    kd1 = 2*PI/sqrt(h11*h11+h12*h12+h13*h13);
    kd2 = 2*PI/sqrt(h21*h21+h22*h22+h23*h23);
    kd3 = 2*PI/sqrt(h31*h31+h32*h32+h33*h33);
    volume = fabs(volume);
    return;
}  /* end init_cell() */


/* the time ratio of evaluating a single term */
/* in the real and reciprocal space summation */
#define TRTF 5.5
/* empirical constants of Ju Li to the Fincham */
/* formulas to achieve the acclaimed accuracy  */
#define RCUT_COEFF 1.2
#define KCUT_COEFF RCUT_COEFF

void init_ewald (int *number_particles,
		 double H[3][3], double *accuracy)
{
    int i,j,k,max_err;
    double xx,yy,zz,r2,maxr;

    NP = *number_particles;
    my_accuracy = *accuracy;
    init_cell(H);
    
    /* Set the parameters alpha and cutoffs based on
       the formula by Fincham CCP5 38, p17 (1993) */ 
    alpha = pow(NP*PI*PI*PI*TRTF/volume/volume,1./6.);
    rcut = sqrt(-log(my_accuracy))/alpha*RCUT_COEFF;
    kcut = 2*alpha*sqrt(-log(my_accuracy))*KCUT_COEFF;

    /* calculate the maximum separation between */
    /* two particles inside the unit cell: must */
    /* be one of the eight corners. */ 
    maxr = 0.;
    for (i=-1; i<=1; i+=2)
	for (j=-1; j<=1; j+=2)
	    for (k=-1; k<=1; k+=2)
	    {
		xx = i*h11 + j*h21 + k*h31;
		yy = i*h12 + j*h22 + k*h32;
		zz = i*h13 + j*h23 + k*h33;
		r2 = xx*xx + yy*yy + zz*zz;
		if (r2>maxr*maxr) maxr=sqrt(r2);
	    }
    
    /* construct a list of important real-space cells */
    maxr += rcut;
    max_xvec = ceil(4.*PI*maxr*maxr*maxr/3./volume*1.2/16)*16;
    nvec = (int *) malloc(3*max_xvec*sizeof(int));
    xvec = (double *) malloc(3*max_xvec*sizeof(double));
    max_n1 = ceil(rcut/rd1);
    max_n2 = ceil(rcut/rd2);
    max_n3 = ceil(rcut/rd3);
    /* first record is the bare unit cell */
    num_xvec = 3;
    xvec[0] = 0.;
    xvec[1] = 0.;
    xvec[2] = 0.;
    for (i=-max_n1; i<=max_n1; i++)
	for (j=-max_n2; j<=max_n2; j++)
	    for (k=-max_n3; k<=max_n3; k++)
		if (!((i==0)&&(j==0)&&(k==0)))
		{
		    xx = i*h11 + j*h21 + k*h31;
		    yy = i*h12 + j*h22 + k*h32;
		    zz = i*h13 + j*h23 + k*h33;
		    r2 = xx*xx + yy*yy + zz*zz;
		    /* only these cells are possible to */
		    /* have an interaction within rcut  */
		    if (r2<maxr*maxr)
		    {
			num_xvec += 3;
			if (num_xvec >= 3*max_xvec)
			{
			    printf ("init_ewald(): max_xvec reached.\n");
			    exit(1);
			}
			nvec[num_xvec-3] = i;
			nvec[num_xvec-2] = j;
			nvec[num_xvec-1] = k;
		    }
		}

    /* construct a list of necessary reciprocal k-points */
    max_gvec = ceil(2.*PI*kcut*kcut*kcut/3./
		    (8*PI*PI*PI/volume)*1.2/16)*16;
    gvec = (int *) malloc(3*max_gvec*sizeof(int));
    max_g1 = ceil(kcut/kd1);
    max_g2 = ceil(kcut/kd2);
    max_g3 = ceil(kcut/kd3);
    /* first record is G=0, which has no inversion image */
    num_gvec = 3;
    gvec[0] = 0;
    gvec[1] = 0;
    gvec[2] = 0;
    /* There are inversion symmetry in energy,  */
    /* force, stress, but not in force constant */
    /* matrix calculations in a general system  */ 
    for (k=0; k<=max_g3; k++)
	for (j=-max_g2; j<=max_g2; j++)
	    for (i=-max_g1; i<=max_g1; i++)
		if ((k>0)||(j>0)||((j==0)&&(i>0)))
		{
		    xx = i*g11 + j*g21 + k*g31;
		    yy = i*g12 + j*g22 + k*g32;
		    zz = i*g13 + j*g23 + k*g33;
		    r2 = xx*xx + yy*yy + zz*zz;
		    if (r2 < kcut*kcut)
		    {
			num_gvec += 3;
			if (num_gvec >= 3*max_gvec)
			{
			    printf ("init_ewald(): max_gvec reached.\n");
			    exit(1);
			}
			gvec[num_gvec-3] = i;
			gvec[num_gvec-2] = j;
			gvec[num_gvec-1] = k;
		    }
		}
    
    /* allocate real space positions and charges */
    x = (double *) malloc(NP*sizeof(double));
    y = (double *) malloc(NP*sizeof(double));
    z = (double *) malloc(NP*sizeof(double));
    q = (double *) malloc(NP*sizeof(double));

#ifdef _TABULATE_ERFC
    /* tabulate the error function and its derivative:     */
    /* because we will do expansion to the 2nd order, the  */
    /* interval should be proportional to (accuracy)^(1/3) */
    derr = pow(*accuracy,1./3.)/2./_TABULATE_ERFC;
    /* default _TABULATE_ERFC = 1. */
    max_err = ceil(alpha*rcut/derr);
    err = (double *)malloc(max_err*sizeof(double));
    cli = (double *)malloc(max_err*sizeof(double));
    for (i=0; i<max_err; i++)
    {
	err[i] = erfc(i*derr);
	cli[i] = 2./sqrt(PI)*exp(-i*i*derr*derr);
    }
#endif

#ifdef _TABULATE_SINE
    /* massive tabulation of sine and cosine values */
    sinn = (double *) malloc(NP*(2*max_g1+1)*(2*max_g2+1)
			     *(2*max_g3+1)*sizeof(double));
    cosn = (double *) malloc(NP*(2*max_g1+1)*(2*max_g2+1)
			     *(2*max_g3+1)*sizeof(double));
#endif

#ifdef _PRINT_EWALD
  printf ("\nAlpha = %f  Rcut = %f  Kcut = %f\n",
	  alpha, rcut, kcut);
  printf ("Max_n1 = %d  Max_n2 = %d  Max_n3 = %d => %d L-points\n",
	  max_n1, max_n2, max_n3, num_xvec/3);
  printf ("Max_g1 = %d  Max_g2 = %d  Max_g3 = %d => %d G-points\n\n",
	  max_g1, max_g2, max_g3, num_gvec/3); 
#endif
  
  return;
} /* end init_ewald() */


void exit_ewald()
{
    free(gvec); free(xvec); free(nvec);
    free(x); free(y); free(z); free(q); 
#ifdef _TABULATE_ERFC 
    free(err); free(cli);
#endif
#ifdef _TABULATE_SINE 
    free(sinn); free(cosn); 
#endif
    return;
} /* end exit_ewald() */


void ewald (charge, s1, s2, s3, H, pote, fx, fy, fz, stress)
/* charge of particles from 1 to NP */
double *charge;
/* reduced coordinates of particles in [-0.5, 0.5] */
double *s1, *s2, *s3;
/* H matrix of the cell, assumed to be passed from fortran */
/* code (column order): in there H(1,1),H(1,2),H(1,3) are  */
/* the xyz components of the first edge, etc. */
double H[3][3]; 
/* Coulomb energy per cell of the system */
double *pote;
/* Coulomb force on particles */
double *fx, *fy, *fz;
/* stress tensor due to Coulomb interactions */
double stress[3][3]; 
{
    double dx,dy,dz,qsum=0.;
    double rx,ry,rz,r2,r,product,xx,margin,factor,gactor;
    double qij,ff;
    double kx,ky,kz,k2,cc,cossum,sinsum,ak;
    double sinx,cosx,siny,cosy,sinz,cosz;
    int i,j,k,l,n;
    int i000,ij00,ijk0,ijkl,jkl;

    /* make it work for Parrinello-Rahman MD */
    init_cell(H);
    alpha = pow(NP*PI*PI*PI*TRTF/volume/volume,1./6.);
    rcut = sqrt(-log(my_accuracy))/alpha*RCUT_COEFF;
    kcut = 2*alpha*sqrt(-log(my_accuracy))*KCUT_COEFF;
    /* the formulas have nice scaling with volume, so */
    /* generated L and G lattices need not be revised */
    /* unless large shape changes happen during MD:   */
    for (n=0; n<num_xvec; n+=3)
    {
	xvec[n]   = nvec[n]*h11 + nvec[n+1]*h21 + nvec[n+2]*h31;
	xvec[n+1] = nvec[n]*h12 + nvec[n+1]*h22 + nvec[n+2]*h32;
	xvec[n+2] = nvec[n]*h13 + nvec[n+1]*h23 + nvec[n+2]*h33;
    }
    /* even the erfc() table has the same range as before */
    
    for (i=0; i<NP; i++)
    {
	fx[i] = 0.;
	fy[i] = 0.;
	fz[i] = 0.;
	while ((s1[i]>0.5)||(s1[i]<-0.5)||
	       (s2[i]>0.5)||(s2[i]<-0.5)||
	       (s3[i]>0.5)||(s3[i]<-0.5))
	{
	    printf("ewald(): reduced coordinates > 0.5.\n");
	    printf("may lead to inaccurate results.\n");
	    /* exit(1); */
	    if (s1[i]<-0.5) s1[i]++;
	    if (s1[i]>0.5)  s1[i]--;
	    if (s2[i]<-0.5) s2[i]++;
	    if (s2[i]>0.5)  s2[i]--;
	    if (s3[i]<-0.5) s3[i]++;
	    if (s3[i]>0.5)  s3[i]--;
	    printf("reduced coordinates modified.\n");
	}
	x[i] = s1[i]*h11 + s2[i]*h21 + s3[i]*h31;
	y[i] = s1[i]*h12 + s2[i]*h22 + s3[i]*h32;
	z[i] = s1[i]*h13 + s2[i]*h23 + s3[i]*h33;
	qsum += charge[i];
    }
    
    if (fabs(qsum/NP) > 0.0001)
	printf("\nwarning from ewald():"\
	       "significant net charge in the system.\n");
    
    for (i=0; i<3; i++) 
	for (j=0; j<3; j++) 
	    stress[i][j] = 0.;
  
    *pote = 0.;
    for (i=0; i<NP; i++)
    {
	q[i] = charge[i] - qsum/NP;
	/* the self-energy */
	*pote -= alpha*q[i]*q[i]/sqrt(PI);
    } 

    /* Do the real space summation */
    for (i=0; i<NP; i++)
	for (j=i; j<NP; j++)
	{
	    dx = x[i] - x[j];
	    dy = y[i] - y[j];
	    dz = z[i] - z[j];
	    qij = q[i] * q[j];
	    for (n=0; n<num_xvec; n+=3)
	    {
		rx = dx + xvec[n];
		ry = dy + xvec[n+1];
		rz = dz + xvec[n+2];
		r2 = rx*rx + ry*ry + rz*rz;
		if ((r2>0.)&&(r2<rcut*rcut))
		{
		    r = sqrt(r2);
		    product = alpha * r;
#ifdef _TABULATE_ERFC
		    k = floor(product/derr);
		    xx = k * derr;
		    margin = product - xx;
		    gactor = err[k] - cli[k] * margin
			* (1-xx*margin);
#else
		    gactor = erfc(product);
#endif
		    /* energy */
		    *pote += (i==j)? gactor*qij/r*0.5:
			gactor*qij/r;
		    
		    /* force */
#ifdef _TABULATE_ERFC
		    margin = product*product - xx*xx;
		    factor = gactor + product * cli[k] *
			(1. - margin * (1 - margin * 0.5));
#else
		    factor = gactor + 2 * product *
			exp(-product*product) / sqrt(PI);
#endif
		    ff = factor * qij / r2 / r;
		    fx[i] += ff * rx;
		    fy[i] += ff * ry;
		    fz[i] += ff * rz; 
		    fx[j] -= ff * rx;
		    fy[j] -= ff * ry;
		    fz[j] -= ff * rz;
		    
		    /* stress */
		    stress[0][0] += (i==j)?ff*rx*rx*0.5:ff*rx*rx;
		    stress[1][1] += (i==j)?ff*ry*ry*0.5:ff*ry*ry;
		    stress[2][2] += (i==j)?ff*rz*rz*0.5:ff*rz*rz;
		    stress[0][1] += (i==j)?ff*rx*ry*0.5:ff*rx*ry;
		    stress[0][2] += (i==j)?ff*rx*rz*0.5:ff*rx*rz;
		    stress[1][2] += (i==j)?ff*ry*rz*0.5:ff*ry*rz;
		}
	    }
	}

    /* Do the reciprocal space summation */

#ifdef _TABULATE_SINE
    /* bootstrap the sine and cosine arrays */
    for (i=0; i<NP; i++)
    {
	sinx = sin(s1[i]*2*PI);
	cosx = cos(s1[i]*2*PI);
	siny = sin(s2[i]*2*PI);
	cosy = cos(s2[i]*2*PI);
	sinz = sin(s3[i]*2*PI);
	cosz = cos(s3[i]*2*PI);

	/* corner */
	i000 = i*(2*max_g1+1)*(2*max_g2+1)*(2*max_g3+1);
	sinn[i000] = sin(2*PI*(-max_g1*s1[i]
			       -max_g2*s2[i]-max_g3*s3[i]));
	cosn[i000] = cos(2*PI*(-max_g1*s1[i]
			       -max_g2*s2[i]-max_g3*s3[i]));

	/* corner -> line */
	for (j=1; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    sinn[ij00] =
		sinn[ij00-(2*max_g2+1)*(2*max_g3+1)]*cosx + 
		cosn[ij00-(2*max_g2+1)*(2*max_g3+1)]*sinx;
	    cosn[ij00] =
		cosn[ij00-(2*max_g2+1)*(2*max_g3+1)]*cosx -
		sinn[ij00-(2*max_g2+1)*(2*max_g3+1)]*sinx;
	}
	
	/* line -> surface */
	for (j=0; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    for (k=1; k<=2*max_g2; k++)
	    {
		ijk0 = ij00 + k*(2*max_g3+1);
		sinn[ijk0] =
		    sinn[ijk0-(2*max_g3+1)]*cosy +
		    cosn[ijk0-(2*max_g3+1)]*siny;
		cosn[ijk0] =
		    cosn[ijk0-(2*max_g3+1)]*cosy -
		    sinn[ijk0-(2*max_g3+1)]*siny;
	    }
	}

	/* surface -> cube */
	for (j=0; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    for (k=0; k<=2*max_g2; k++)
	    {
		ijk0 = ij00 + k*(2*max_g3+1);
		for (l=1; l<=2*max_g3; l++)
		{
		    ijkl = ijk0 + l;
		    sinn[ijkl] =
			sinn[ijkl-1]*cosz +
			cosn[ijkl-1]*sinz;
		    cosn[ijkl] =
			cosn[ijkl-1]*cosz -
			sinn[ijkl-1]*sinz;
		}
	    }
	}
    }  
#endif
    
    /* in the summation, omit G = 0 */
    for (n=3; n<num_gvec; n+=3)
    {
	cossum = 0.;
	sinsum = 0.;
	
#ifdef _TABULATE_SINE
	jkl = ((gvec[n]+max_g1)*(2*max_g2+1)+ 
	        gvec[n+1]+max_g2)*(2*max_g3+1) +
	        gvec[n+2]+max_g3;
#endif

	for (i=0; i<NP; i++)
	{
#ifdef _TABULATE_SINE
	    ijkl = i*(2*max_g1+1)*(2*max_g2+1)
		*(2*max_g3+1) + jkl;
	    cossum += q[i]*cosn[ijkl];
	    sinsum += q[i]*sinn[ijkl];
#else
	    cossum += q[i]*cos(2*PI*(gvec[n]*s1[i] + 
				     gvec[n+1]*s2[i] +
				     gvec[n+2]*s3[i]));
	    sinsum += q[i]*sin(2*PI*(gvec[n]*s1[i] +
				     gvec[n+1]*s2[i] +
				     gvec[n+2]*s3[i]));
#endif
	}

	kx = gvec[n]*g11+gvec[n+1]*g21+gvec[n+2]*g31;
	ky = gvec[n]*g12+gvec[n+1]*g22+gvec[n+2]*g32;
	kz = gvec[n]*g13+gvec[n+1]*g23+gvec[n+2]*g33;
	k2 = kx*kx + ky*ky + kz*kz;
	/* with inversion symmetry, each G represents two */
	ak = 4.*PI*exp(-k2/alpha/alpha/4.)/k2/volume;

	/* energy */
	cc = ak*(cossum*cossum+sinsum*sinsum);
	*pote += cc;

	/* force */
	for (i=0; i<NP; i++)
	{
	    
#ifdef _TABULATE_SINE             
	    ijkl = i*(2*max_g1+1)*(2*max_g2+1)
		*(2*max_g3+1) + jkl;
	    cc = 2.*ak*q[i]*(cossum*sinn[ijkl]-
			     sinsum*cosn[ijkl]);
#else
	    cc = 2*ak*q[i]*(cossum*sin(2*PI*(gvec[n]*s1[i] +
					     gvec[n+1]*s2[i] +
					     gvec[n+2]*s3[i])) -
			    sinsum*cos(2*PI*(gvec[n]*s1[i] +
					     gvec[n+1]*s2[i] +
					     gvec[n+2]*s3[i])));
#endif
	    fx[i] += kx*cc;
	    fy[i] += ky*cc;
	    fz[i] += kz*cc;
	}
    }

    return;    
} 
