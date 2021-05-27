#!/bin/bash 


for i in $(seq 1 8)
do
gcc -fopenmp prog_in.c -o in_$i.x -lm
export OMP_NUM_THREADS=$i
( time ./in_$i.x)2> time_$i.dat
done 
