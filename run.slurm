#!/bin/bash
#SBATCH -n 96
#SBATCH -q debug
#SBATCH -t 01:00:00
srun -n 48 -c2 --cpu_bind=cores ./sa.x control.PT > simulated.out
