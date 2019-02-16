#!/usr/bin/env python
import random
import sys
import math
import numpy as np
def nvt(t,tick,file,step):
	file.write("fix "+str(tick)+" all"+" nvt"+" temp"+" "+str(t)+" "+str(t)+" "+"1.0"+"\n");
	file.write("run"+" "+str(step)+"\n");
	file.write("unfix"+" "+str(tick)+"\n");
	file.write("\n");
def npt(t,tick,file,step):
	file.write("fix "+str(tick)+" all"+" npt"+" temp"+" "+str(t)+" "+str(t)+" "+"1.0"+" "+"aniso"+" "+"1.01325 "+"1.01325"+" 5.0\n");
	file.write("run "+str(step)+"\n");
	file.write("unfix "+str(tick)+"\n");
	file.write("\n");
def dump(t,tick,file,step):
	file.write("fix "+str(tick)+" all"+" npt"+" temp"+" "+str(t)+" "+str(t)+" "+"1.0"+" "+"aniso"+" "+"1.01325 "+"1.01325"+" 5.0\n");
	file.write("dump "+str(tick)+" all custom 200 "+"dump"+str(t)+".xyz"+" "+"x y z\n");
	file.write("dump_modify "+str(tick)+" sort"+" id\n");
	file.write("run "+str(step)+"\n");
	file.write("undump "+str(tick)+"\n");
	file.write("unfix "+str(tick)+"\n");
	file.write("\n");
infile="test.lammpsinput";
outfile="new.BZO";
fin=open(infile,"r");
fout=open(outfile,"w");
for i in range(55):
	temp=fin.readline();
	fout.write(temp);
tick=1;
listone=range(10,100,10);
listtwo=range(120,320,40);
listall=listone+listtwo;
count=0;
equi=40000;
dum=600000;
for i in listall:
    count=count+2*equi+dum;
    nvt(i,tick,fout,equi);
    tick=tick+1;
    npt(i,tick,fout,equi);
    tick=tick+1;
    dump(i,tick,fout,dum);
print count;
