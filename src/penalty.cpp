#include "penalty.h"
#include <iomanip>
#include <new>
#include "atom.h"
#include "readpara.h"
#include "readion.h"
#include <math.h>
#include <mpi.h>
#include "sa.h"
/*Overall comment by Zhenbang: to ensure that this code can work, say we have 
4 species, then site should be 4, charge should be 4, xp charge part should be 3
*/

/*
 * the main AIM of the two function is to map the list variable into bvv model parameters and charge parameters. 
 *
 * /
/*tell you which index you are going to map tick to (i,j) */
void indexbvvmap(int** bvvmatrixmap,int tick,int& i,int& j){
	size_t count=0;
	for(size_t m=0;m<control::pair_num;m++)
		for(size_t n=0;n<12;n++){
			if(bvvmatrixmap[m][n]==1){
				if(count==tick){
					i=m;
					j=n;
				}
				count++;
			}
		}
}
/*this tell you which index you are going to map tick to (i) in chargemap
*now we are dealing with only asite, bsite and osite charge change.
tick is the index of xp in the map function
*/
void indexchargemap(int* chargemap,int tick,int& i){
	int count=tick+1;
	int remain=count-control::paracount_bvv;/*the remaining parameters to be optimized, Assumed to be in the charge map*/
	size_t sum=0;
	for(size_t m=0;m<species::spe.size();m++){
		if(chargemap[m]==1){
			sum=sum+1;
			if(sum==remain){
				i=m;
			}
		}
	}
}
void mapjiahao(double* xp){
	double sum=0.0;
	for(size_t i=0;i<control::paracount_bvv;i++){
		control::bvvmatrix[control::mapXpTickToBvvTick[i][1]][control::mapXpTickToBvvTick[i][2]]=xp[i];
	}
	for(size_t i=control::paracount_bvv;i<control::independentcharge+control::paracount_bvv;i++){
		control::chargexp[control::mapXpTickToChargeTick[i-control::paracount_bvv][1]]=xp[i];
	}
  size_t size=species::spe.size();
  for(size_t i=0;i<size;i++){
    sum=0.0;
    for(size_t j=0;j< control::independentcharge;j++){
      sum=sum+control::chargemap[i][j]*control::chargexp[j];
    }
    control::charge[i]=sum+control::chargeres[i];
  }
}
double PenaltyFunc(double* xp, box* system,int numberone, int index,int databasetick){
    int indexRef = index;
    mapjiahao(xp);
    double penalty = 0.0;
    box* ionall = system; 
    int number = numberone;
		int world_rank,mpi_size;
		MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
		MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
		int box_ave=floor((number+0.00001)/mpi_size);
		int remain=number%mpi_size;
		int box_size_local;
		double ref_energy[2]={0.0,0.0};
		if(world_rank<remain){
			box_size_local=box_ave+1;
		}
		else{
			box_size_local=box_ave;
		}
		int ref_proc=indexRef%mpi_size;
		double* mdenergy=new double [number];
		double* dftenergy=new double [number];
		double* diffenergy=new double [number];
		for(size_t i=0;i<number;i++){
			mdenergy[i]=0.0;
			dftenergy[i]=0.0;
			diffenergy[i]=0.0;
		}
		int ref_id_proc=floor((indexRef+0.000)/mpi_size);
    for (size_t i=0; i<box_size_local; i++){
        ionall[i].mdenergy = 0;
        ionall[i].updatebvparameter(control::bvvmatrix);
        ionall[i].computeAll();
    }
		if(ref_proc==world_rank){
			ref_energy[0]=ionall[ref_id_proc].mdenergy;			
			ref_energy[1]=ionall[ref_id_proc].dftenergy;
//			std::cout<<"the reference energy is: "<<ref_energy[0]<<" "<<std::setprecision(15)<<ref_energy[1]<<std::endl;
		}
		MPI_Bcast(&ref_energy,2,MPI_DOUBLE,ref_proc,MPI_COMM_WORLD);
    /*Calculate the penalty*/
    double PenaltyE = 0;
    double PenaltyF = 0;
    for (size_t i=0; i<box_size_local; i++){
        PenaltyE += fabs((ionall[i].mdenergy-ref_energy[0]) - (ionall[i].dftenergy-ref_energy[1]))*ionall[i].weight;
				diffenergy[i*mpi_size+world_rank]=fabs((ionall[i].mdenergy-ref_energy[0]) - (ionall[i].dftenergy-ref_energy[1]));
				dftenergy[i*mpi_size+world_rank]=ionall[i].dftenergy-ref_energy[1];
				mdenergy[i*mpi_size+world_rank]=ionall[i].mdenergy-ref_energy[0];
				for (size_t j=0; j<ionall[i].size; j++){
            for (size_t k=0; k<3; k++){
                PenaltyF += fabs((ionall[i].allatom[j].force[k]-ionall[i].allatom[j].dftforce[k]))*ionall[i].weight; 
						}
        }
    }
    PenaltyE = PenaltyE/number*saconst::sa_eweight;
    PenaltyF = PenaltyF/(number*3*ionall[0].size)*saconst::sa_fweight;
		penalty = PenaltyE + PenaltyF;
		double Penaltyall=0.0;
		MPI_Reduce(&penalty,&Penaltyall,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(mdenergy,control::mdenergy[databasetick],number,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(dftenergy,control::dftenergy[databasetick],number,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Reduce(diffenergy,control::diffenergy[databasetick],number,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		delete [] mdenergy;
		delete [] dftenergy;
		delete [] diffenergy;
    return Penaltyall;
}
