#!/usr/bin/env python2
import numpy as np
import sys
crystal_coord=np.loadtxt("posit3");
crystal_length=[8.0276629746782101,8.0406596742829339,8.0199791338267016];
atom=(crystal_coord.shape)[0];
dim=(crystal_coord.shape)[1];
for i in range(atom):
	for j in range(dim):
		crystal_coord[i][j]=crystal_coord[i][j]*crystal_length[j];
np.savetxt(sys.stdout, crystal_coord, '%15.10f')
