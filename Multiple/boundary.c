#include "boundary.h"
#include <stdio.h>

//grid is parallelised in the x direction

void boundarypsi(double **psi, int m, int n,
		 int b, int h, int w, MPI_Comm comm)
{

  int size, rank,i,j;
  int istart, istop;

  MPI_Comm_size(comm,&size);
  MPI_Comm_rank(comm,&rank);

  istart = m*rank + 1;
  istop = istart + m -1;

  //BCs on bottom edge

  for (i=b+1;i<=b+w-1;i++)
    {
      if (i >= istart && i <= istop)
	{
	  psi[i-istart+1][0] = (float)(i-b);
        }
    }

  for (i=b+w;i<m*size+1;i++)
    {
      if (i >= istart && i <= istop)
	{
	  psi[i-istart+1][0] = (float)(w);
	}
    }

  //BCS on RHS

  if (rank == size-1)
    {
      for (j=1; j <= h; j++)
	{
	  psi[m+1][j] = (float) w;
	}

      for (j=h+1;j<=h+w-1; j++)
	{
	  psi[m+1][j]=(float)(w-j+h);
	}
    }
}

void boundaryzet(double **zet, double **psi, int m, int n, MPI_Comm comm)
{
  int size, rank, i,j;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  //set top/bottom BCs:

  for (i=1;i<m+1;i++)
    {
      zet[i][0] = 2.0*(psi[i][1]-psi[i][0]);
      zet[i][n+1] = 2.0*(psi[i][n]-psi[i][n+1]);
    }

  //set left BC:
  if (rank ==0)
    {
      for (j=1;j<n+1;j++)
	{
	  zet[0][j] = 2.0*(psi[1][j]-psi[0][j]);
	}
    }

  //set left BCs

  if (rank ==size-1)
    {
      for (j=1;j<n+1;j++)
	{
	  zet[m+1][j] = 2.0*(psi[m][j]-psi[m+1][j]);
	}
    }
}

void haloswap_thread(double **x, int m, int n, MPI_Comm comm, int myid)
{
  int rank, uprank, dnrank, size;
  //int tag=1;

  MPI_Status status;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);


  //no need to halo swap if serial:

  if (size > 1)
    {
      // determine uprank and downrank
      if(rank == size-1)
	{
	  uprank = MPI_PROC_NULL;
        }
      else
	{
	  uprank = rank+1;
        }

      if (rank == 0)
	{
	  dnrank = MPI_PROC_NULL;
        }
      else
	{
	  dnrank = rank-1;
        }

      //send right boundaries and receive left ones

      MPI_Sendrecv(&x[m][1+n*myid],n,MPI_DOUBLE,uprank,myid,
		   &x[0][1+n*myid],n,MPI_DOUBLE,dnrank,myid,
		   comm,&status);

      //send left boundary and receive right

      MPI_Sendrecv(&x[1][1+n*myid],  n,MPI_DOUBLE,dnrank,myid,
		   &x[m+1][1+n*myid],n,MPI_DOUBLE,uprank,myid,
		   comm,&status);

      //printf("Swapped everything :) \n");
    }
}

void haloswap(double **x, int m, int n, MPI_Comm comm)
{
  int rank, uprank, dnrank, size;
  int tag=1;

  MPI_Status status;

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);


  //no need to halo swap if serial:

  if (size > 1)
    {
      // determine uprank and downrank
      if(rank == size-1)
        {
          uprank = MPI_PROC_NULL;
        }
      else
        {
          uprank = rank+1;
        }

      if (rank == 0)
        {
          dnrank = MPI_PROC_NULL;
        }
      else
        {
          dnrank = rank-1;
        }

      //send right boundaries and receive left ones

      MPI_Sendrecv(x[m],n,MPI_DOUBLE,uprank,tag,
                   x[0],n,MPI_DOUBLE,dnrank,tag,
                   comm,&status);

      //send left boundary and receive right

      MPI_Sendrecv(x[1],n,MPI_DOUBLE,dnrank,tag,
                   x[m+1],n,MPI_DOUBLE,uprank,tag,
                   comm,&status);

      //printf("Swapped everything :) \n");
    }
}

