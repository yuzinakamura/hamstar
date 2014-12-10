#!/bin/bash

#set eps, sweep time
#4 16 32
for NUM_PROC in 1 8
do
  rm stats_${NUM_PROC}.csv
  for SEED in 1 2 3 4
  do
    mpic++ -std=c++0x -o SMHA SMHA.cpp
    mpirun -n ${NUM_PROC} -hostfile host_file ./SMHA -seed ${SEED}
  done
done