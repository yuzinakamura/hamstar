#!/bin/bash

#set eps, sweep time
#4 16 32
mpic++ -std=c++0x -o SMHA SMHA.cpp

for NUM_PROC in 1 8
do
  rm stats_${NUM_PROC}.csv
  for SEED in 1 2 3 4
  do
    mpirun -n ${NUM_PROC} -hostfile host_file ./SMHA -seed ${SEED} -filename stats_${NUM_PROC}
  done
done