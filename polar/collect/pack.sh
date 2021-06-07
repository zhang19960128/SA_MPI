#!/bin/bash
path=$(pwd);
all.dat
for i in $(seq 1 5)
do
    g++ database.cpp -o col
    ./col relax.out$i dat 0
    cat dat >>all.dat
done
