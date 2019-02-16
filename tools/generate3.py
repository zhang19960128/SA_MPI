#!/usr/bin/env python
import random
import sys
import math
import numpy as np
'''
angle type
here 0 BBBB
     1 BBBC
		 2 BBCC
		 3 BCBC
		 4 CCCB
		 5 CCCC
'''
def angle(inp):
	if(np.sum(inp))==4:
		return 5+1;
	elif(np.sum(inp)==3):
		return 4+1;
	elif(np.sum(inp)==1):
		return 1+1;
	elif(np.sum(inp)==0):
		return 0+1;
	elif(inp[0]==inp[2]):
		return 3+1;
	else:
		return 2+1;
def changeindex(index,cell):
	'''change the index to 3D array[index_x,index_y,index_z]'''
	re=np.zeros(3);
	re[2]=math.floor(index/cell/cell);
	re[1]=math.floor((index-re[2]*cell*cell)/cell);
	re[0]=index-re[2]*cell*cell-re[1]*cell;
	return re;
def backindex(indexlist,cell):
	'''change the three dimensional index into a one dimentional one'''
	return indexlist[0]%cell+(indexlist[1]%cell)*cell+(indexlist[2]%cell)*cell*cell;
def typein(index,Asite):
	a=[];
	for i in index:
		a.append(Asite[int(i)]);
	return a;
fraction=1.0;
cell=2;'''here are the period on x y z cell length'''
ba_valen=sys.argv[1];
ca_valen=sys.argv[1];
zr_valen=sys.argv[2];
ox_valen=sys.argv[3];
calist=random.sample(range(cell*cell*cell),int(cell*cell*cell*fraction));
calist=sorted(calist);
ca=open("cadata.txt","w");
for i in calist:
	ca.write(str(i)+"\n");
ca.close();
data=open("mixdata.BTO","w");
data.write("#LAMMPS 100 water\n");
data.write("\n");
data.write(str(5*cell*cell*cell)+" atoms\n");
data.write(str(3*cell*cell*cell)+" angles\n");
data.write("3 atom types\n");
data.write("6 angle types\n");
data.write("\n");
cell_length=8.386494280/2;
data.write("0.0 "+str(cell_length*cell)+" "+"xlo xhi\n");
data.write("0.0 "+str(cell_length*cell)+" "+"ylo yhi\n");
data.write("0.0 "+str(cell_length*cell)+" "+"zlo zhi\n");
data.write("\n");
data.write("\n");
data.write("Masses\n");
data.write("\n");
data.write("1 137.327  #Ba\n");
data.write("2 91.224 #Zr\n");
data.write("3 15.9993  #o\n");
data.write("\n");
data.write("\n")
'''
angle type
here 0 BBBB 1.37953
     1 BBBC 1.11925
		 2 BBCC 0.869524
		 3 BCBC 0.854168
		 4 CCCB 0.599274
		 5 CCCC 0.33389
'''
data.write("Angle Coeffs\n");
data.write("\n");
data.write("1 "+str(0.000)+" 180.00"+"\n");
data.write("2 "+str(0.000)+" 180.00"+"\n");
data.write("3 "+str(0.000)+" 180.00"+"\n");
data.write("4 "+str(0.000)+" 180.00"+"\n");
data.write("5 "+str(0.000)+" 180.00"+"\n");
data.write("6 "+str(0.000)+" 180.00"+"\n");
data.write("\n");
data.write("Atoms\n");
data.write("\n");
Asite=np.zeros(cell*cell*cell);
for i in calist:
	Asite[i]=1;
for i in range(cell*cell*cell):
	if(Asite[i]==1):
		index=changeindex(i,cell);
		data.write(str(i+1)+" "+"1"+" "+"1"+" "+str(ca_valen)+" "+str(cell_length*index[0])+" "+str(index[1]*cell_length)+" "+str(index[2]*cell_length)+"\n");
	else:
		index=changeindex(i,cell);
		data.write(str(i+1)+" "+"1"+" "+"1"+" "+str(ca_valen)+" "+str(cell_length*index[0])+" "+str(index[1]*cell_length)+" "+str(index[2]*cell_length)+"\n");
for i in range(cell*cell*cell):
	index=changeindex(i,cell);
	data.write(str(i+cell*cell*cell*1+1)+" "+"1"+" "+"2"+" "+str(zr_valen)+" "+str(cell_length*index[0]+cell_length/2)+"  "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2]+cell_length/2)+"\n");
for i in range(cell*cell*cell):
	'''bottom to up'''
	index=changeindex(i,cell);
	data.write(str(i+2*cell*cell*cell+1)+" "+"1"+" "+"3"+" "+str(ox_valen)+" "+str(cell_length*index[0]+cell_length/2)+" "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2])+"\n");
for i in range(cell*cell*cell):
	'''front to back'''
	index=changeindex(i,cell);
	data.write(str(i+3*cell*cell*cell+1)+" "+"1"+" "+"3"+" "+str(ox_valen)+" "+str(cell_length*index[0]+cell_length/2)+" "+str(cell_length*index[1])+" "+str(cell_length*index[2]+cell_length/2)+"\n");
for i in range(cell*cell*cell):
	'''left to right'''
	index=changeindex(i,cell);
	data.write(str(i+4*cell*cell*cell+1)+" "+"1"+" "+"3"+" "+str(ox_valen)+" "+str(cell_length*index[0])+" "+str(cell_length*index[1]+cell_length/2)+" "+str(cell_length*index[2]+cell_length/2)+"\n");
data.write("\n");
data.write("Angles\n");
data.write("\n");
four=np.zeros(4);
for i in range(cell*cell*cell):
	index=changeindex(i,cell);
	atom1=backindex([index[0],index[1],index[2]],cell);
	atom2=backindex([index[0],index[1],(index[2]+1)%cell],cell);
	atom3=backindex([index[0],index[1],(index[2]+2)%cell],cell);
	four[0]=backindex([index[0],index[1],(index[2]+1)],cell);
	four[1]=backindex([index[0]+1,index[1],index[2]+1],cell);
	four[2]=backindex([index[0]+1,index[1]+1,index[2]+1],cell);
	four[3]=backindex([index[0],index[1]+1,index[2]+1],cell);
	atom1=int(atom1);
	atom2=int(atom2);
	atom3=int(atom3);
	data.write(str(i+1)+" "+str(angle(typein(four,Asite)))+" "+str(atom1+1+2*cell*cell*cell)+" "+str(atom2+1+2*cell*cell*cell)+" "+str(atom3+1+2*cell*cell*cell)+"\n");
for i in range(cell*cell*cell):
	index=changeindex(i,cell);
	atom1=backindex([index[0],index[1],index[2]],cell);
	atom2=backindex([index[0],(index[1]+1)%cell,index[2]],cell);
	atom3=backindex([index[0],(index[1]+2)%cell,index[2]],cell);
	atom1=int(atom1);
	atom2=int(atom2);
	atom3=int(atom3);
	four[0]=backindex([index[0],index[1]+1,index[2]],cell);
	four[1]=backindex([index[0]+1,index[1]+1,index[2]],cell);
	four[2]=backindex([index[1]+1,index[1]+1,index[2]+1],cell);
	four[3]=backindex([index[1],index[1]+1,index[2]+1],cell);
	data.write(str(i+cell*cell*cell+1)+" "+str(angle(typein(four,Asite)))+" "+str(atom1+3*cell*cell*cell+1)+" "+str(atom2+3*cell*cell*cell+1)+" "+str(atom3+3*cell*cell*cell+1)+"\n");
for i in range(cell*cell*cell):
	index=changeindex(i,cell);
	atom1=backindex([index[0],index[1],index[2]],cell);
	atom2=backindex([(index[0]+1)%cell,index[1],index[2]],cell);
	atom3=backindex([(index[0]+2)%cell,index[1],index[2]],cell);
	atom1=int(atom1);
	atom2=int(atom2);
	atom3=int(atom3);
	four[0]=backindex([index[0]+1,index[1],index[2]],cell);
	four[1]=backindex([index[0]+1,index[1]+1,index[2]],cell);
	four[2]=backindex([index[0]+1,index[1]+1,index[2]+1],cell);
	four[3]=backindex([index[0]+1,index[1],index[2]+1],cell);
	data.write(str(i+2*cell*cell*cell+1)+" "+str(angle(typein(four,Asite)))+" "+str(atom1+4*cell*cell*cell+1)+" "+str(atom2+4*cell*cell*cell+1)+" "+str(atom3+4*cell*cell*cell+1)+"\n");
