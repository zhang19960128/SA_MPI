#include "penalty.h"
#include <iomanip>
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

void map(double* xp){
	int m,n;
 	for(int i=0;i<control::paracount_bvv;i++){
		indexbvvmap(control::bvvmatrixmap,i,m,n);
		control::bvvmatrix[m][n]=xp[i];
	}

    double accCharge = 0;
    double posAccCharge = 0;
    /*Zhenbang: ensure the same site has the same charge*/ 
    double* siteID = new double[species::spe.size()];
    for (size_t i=0; i<species::spe.size(); i++){
        siteID[i] = -1;
    }
    int count = 0;
    /*store the sites that have appeared. After the loop, count should be either 2 or 3*/
    for(int i=control::paracount_bvv;i<control::paracount_bvv+control::paracount_charge;i++){
				indexchargemap(control::chargemap,i,m);
        control::charge[m]=xp[i];
        int count2=0;
        int sameSite = -1;
        for (int j=0; j<=count; j++){ 
            if (species::site[m] != siteID[j]){
                count2 += 1;
             }
             else{
                sameSite = siteID[j];
             }
        }

        if (count2 == 0){
            siteID[count] = species::site[m];
            count += 1;
        }
        else{
            for (int j=0; j<species::site.size();j++){
                if (sameSite == species::site[j] && m != j){
                    control::charge[m] = control::charge[j];
                }
            }
        }

        if (species::site[m] == 0 || species::site[m] == 1){
            if (count2 == 0){
                accCharge += control::charge[m];//avoid overcounting
                posAccCharge += control::charge[m];
            }
        }
        else{
            accCharge += control::charge[m]*3;
        }
    }
	/*perserve charge conservation law! Zhenbang*/
    for (int i=0; i<species::spe.size(); i++){
        if (control::chargemap[i] == 0 && count == 2){
            if (species::site[i] == 0 || species::site[i] == 1){
                control::charge[i] = (0-accCharge);
            }
            else{
                control::charge[i] = (0-accCharge)/3;
            }
        }
        /*All the site has already been updated. Then take O site as the dependent variable*/
        else if (control::chargemap[i] == 0 && count == 3){
            int thisSite = species::site[i];
            for (int j=0; j<species::spe.size(); j++){
                if (thisSite == species::site[j] && j!= i) control::charge[i] = control::charge[j];
            }
            for (int j=0; j<species::spe.size(); j++){
                if (species::site[j] == 2) control::charge[j] = (0-posAccCharge)/3; 
            }
        }
    }
    delete[] siteID;
}

//Zhenbang
void mapToXp(double* xp){
    int count = 0;
    for (size_t m=0; m<control::pair_num; m++){
        for (size_t n=0; n<12; n++){
            if (control::bvvmatrixmap[m][n] == 0) continue;
            xp[count] = control::bvvmatrix[m][n];
            count += 1;
        }
    }

    for (size_t i=0; i<species::spe.size(); i++){
        if (control::chargemap[i] == 1){
            xp[count] = control::charge[i];//the number of this charge array is Nspecies-1 
            count += 1;
        }
    }
}
/********************************************end by charge************************************/
//Zhenbang
/*
numberone: give how many structures are in the box
index: give the tick of minimum energy in this database.
*/
void mapjiahao(double* xp){
	double sum=0.0;
	int size=control::para_site_charge.size();
	for(size_t i=0;i<control::paracount_bvv;i++){
		control::bvvmatrix[control::mapXpTickToBvvTick[i][1]][control::mapXpTickToBvvTick[i][2]]=xp[i];
	}
	for(size_t i=control::paracount_bvv;i<control::paracount_charge+control::paracount_bvv;i++){
		control::para_site_charge[control::mapXpTickToChargeTick[i-control::paracount_bvv][1]]=xp[i];
	}
	if(control::neutral==1){
		/*force charge neutral*/
	//	control::para_site_charge[control::lastchargetick]=
			for(size_t i=0;i<size;i++){
				if(i==control::lastchargetick){
					continue;
				}
				sum=sum+control::chemical_formula[i]*control::para_site_charge[i];
			}
			control::para_site_charge[control::lastchargetick]=(0-sum)/control::chemical_formula[control::lastchargetick];
			size=species::spe.size();
			for(size_t i=0;i<size;i++){
				control::charge[i]=control::para_site_charge[species::site[i]];
			}
	}
	else{
	};

}
double PenaltyFunc(double* xp, box* system,int numberone, int index){
    int indexRef = index;
    mapjiahao(xp);
    double penalty = 0.0;
    box* ionall = system; 
    int number = numberone;
    double* totalEnergy = new double[number];
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
	//	if(world_rank==0){
	//		std::cout<<"The total Penalty is: "<<Penaltyall<<std::endl;
	//	}
    return Penaltyall;
}
