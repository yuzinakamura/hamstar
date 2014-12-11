Once a node is booted. Run the following command to install OpenMPI:

sudo apt-get install openmpi-bin openmpi-common openssh-client openssh-server libopenmpi1.3 libopenmpi-dbg libopenmpi-dev

The host_file specifies the machines and the number of threads available on each machine.

To compile the source code, run the following command:

mpic++ SMHA.cpp -std=c++0x -o smha -O3

To run it with 32 processes (when only 16 are available on a single machine), do the following:

mpirun -n 32 -hostfile host_file ./smha

It will fill in processes in the order of host_file, filling up each node until the maximum is reached.

5x5, 1 core, 1 seed


Output Format:
Basically stats, then whatever's new, and whatever this particular run is
stats_cores_4.csv