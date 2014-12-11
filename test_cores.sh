#!/bin/bash

mpic++ -std=c++0x -o SMHA SMHA.cpp -O3

for NUM_PROC in 2 4
do
  rm stats_cores10k_${NUM_PROC}.csv
  for SEED in 1 2 3 4 5 6 7 8 9 10
  do
    mpirun -n ${NUM_PROC} -hostfile host_file16 ./SMHA -s ${SEED} -f stats_cores10k_${NUM_PROC}
  done
done