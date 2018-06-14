# Hybrid_MPI-OpenMP
A Jacobi algorithm solution applied by hybrid programming with MPI and OpenMP

cfd.pbs is for Archer and t2fd.pbs is for Cirrus.
The number of threads is written in the cfd.c and it is initialised with 4. We can change the nthreads in cfd.c to change the number of threads. At last, it will be revised and added to cfd.pbs.
