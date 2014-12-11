#!/bin/bash

g++ -std=c++0x -o lame_SMHA lame_SMHA.cpp -O3


for NUM_PROC in 32
do
  rm stats_lame_${NUM_PROC}.csv
  for SEED in 1 2 3 4 5 6 7 8 9 10
  do
    ./lame_SMHA -n ${NUM_PROC} -s ${SEED} -f stats_lame_${NUM_PROC}
  done
done