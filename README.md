This code is built by the people in [Rappe Group](http://web.sas.upenn.edu/rappegroup/)
Contributors are **Jiahao Zhang**, **ZhengBang Dai** under the supervision of **Prof. Rappe**.  
If you encounter any problems, please don't contact us.  
1. Simulated Annealing Code for Bond-Valence parameters generation.
2. Fully paralleled with MPI code.
3. Every proc are govern with it's own database region.
4. Capable of generating the same parameter but different database configurations.
5. For DOD thunder machines:  
			module load intel-mpi  
			Change Makefile **CXX=##** to **CXX=mpiicc**
			module rm gnu-compiler  
			make -j 8  
6. For DOE machine:  
			Change Makefile **Cxx=##** to **Cxx=CC**

