import numpy as np
a=[8.0411057168277760,0.0,0.0];
b=[0.5,8.0411057168277760,0.0];
c=[0.0,0.0,8.0411057168277760];
crystal=np.loadtxt("posit1");
natom=40;
for i in range(natom):
  re=np.add(np.add(np.multiply(crystal[i][0],a),np.multiply(crystal[i][1],b)),np.multiply(crystal[i][2],c));
  print(str(re[0])+" "+str(re[1])+" "+str(re[2]))
