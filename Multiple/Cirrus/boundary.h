#include <mpi.h>

void boundarypsi(double **psi, int m, int n, int b, int h, int w,
		 MPI_Comm comm);

void boundaryzet(double **zet, double **psi, int m, int n, MPI_Comm comm);

void haloswap_thread(double **x, int m, int n, MPI_Comm comm, int myid);

void haloswap(double **x, int m, int n, MPI_Comm comm);
