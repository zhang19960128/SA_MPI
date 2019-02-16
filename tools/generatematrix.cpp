#include <iostream>
#include <new>
int main(){
	int elements=3;//specify how many elements are in the simulation
	int pair=elements*(elements+1)/2;//how many pairs iteractions are in the simulation
	int newelements=2;//which is the new elements. the parameter don't know, the tick go from 0 1 2 3 4 5 6 7
	int** matrix=new int* [pair];
	for(size_t i=0;i<pair;i++){
		matrix[i]=new int [12];
	}
	size_t count=0;
	for(size_t i=1;i<=elements;i++)
		for(size_t j=i;j<=elements;j++){
			/*specify the bond valence itself*/
			if(i==newelements && j==newelements){
				matrix[count][0]=0;
				matrix[count][1]=0;
				matrix[count][2]=1;
				matrix[count][3]=0;
				matrix[count][4]=0;
				matrix[count][5]=0;
				matrix[count][6]=0;
				matrix[count][7]=1;
				matrix[count][8]=0;
				matrix[count][9]=0;
				matrix[count][10]=1;
				matrix[count][11]=1;
			}
      /*specify the bond valence between this guy and the negative cations, i is the element, j is the anions*/
			else if( i==newelements && j==elements){
				matrix[count][0]=1;
				matrix[count][1]=1;
				matrix[count][2]=0;
				matrix[count][3]=0;
				matrix[count][4]=0;
				matrix[count][5]=0;
				matrix[count][6]=0;
				matrix[count][7]=1;
				matrix[count][8]=0;
				matrix[count][9]=0;
				matrix[count][10]=0;
				matrix[count][11]=0;
			}
			/*the only thing that change with new cations and other cations in the lj force function*/
			else if( j!=elements && i!=elements && (i!=j) && (i==newelements || j==newelements)){
				matrix[count][0]=0;
				matrix[count][1]=0;
				matrix[count][2]=0;
				matrix[count][3]=0;
				matrix[count][4]=0;
				matrix[count][5]=0;
				matrix[count][6]=0;
				matrix[count][7]=1;
				matrix[count][8]=0;
				matrix[count][9]=0;
				matrix[count][10]=0;
				matrix[count][11]=0;
			}
			else{
			for(size_t k=0;k<12;k++){
				matrix[count][k]=0;
			}
			}
			count++;
		}
	for(size_t i=0;i<pair;i++){
		for(size_t j=0;j<12;j++){
			std::cout<<matrix[i][j]<<"\t";
		}
		std::cout<<std::endl;
	}
}
