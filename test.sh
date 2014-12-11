#!/bin/bash

#set eps, sweep time
#4 16 32
mpic++ -std=c++0x -O3 -o SMHA SMHA.cpp
g++ -std=c++0x -O3 -o lame_SMHA lame_SMHA.cpp


for SEED in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
  rm stats_cool_${NUM_PROC}.csv
  rm stats_lame_${NUM_PROC}.csv
  for NUM_PROC in 1 2 4 8 16
  do
    mpirun -n ${NUM_PROC} -hostfile host_file ./SMHA -s ${SEED} -f stats_cool_${NUM_PROC}
    ./lame_SMHA -n ${NUM_PROC} -s ${SEED} -f stats_lame_${NUM_PROC}
  done
done