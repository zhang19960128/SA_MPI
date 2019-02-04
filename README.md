1. Simulated Annealing Code for Bond-Valence parameters generation.
2. Fully paralleled with MPI code.
3. Every proc are govern with it's own database region.
4. Capable of generating the same parameter but different database configurations.
5. For DOD thunder machines:  
			module load intel-mpi  
			module rm gnu-compiler  
			make -j 8  
6. For DOE machine:  
			module load impi  
			make -j 8  
