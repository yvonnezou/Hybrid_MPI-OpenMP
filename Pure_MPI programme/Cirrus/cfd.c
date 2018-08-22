/**
 *
 * This the pure-MPI programme on Cirrus. We add a timer in it and remove some useless code to keep the output clean.
 * Also some code in cfdio.c has been removed for the same reason.
 * The difference between this programme and the one on ARCHER is the size of time_iter, time_comp and time_comm.
 * The pbs file is also different because of different wrappers and MPI libraries.
 * We add two shell scripts run.sh and runsplit.sh.
 *
**/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

#include "arraymalloc.h"
#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double localerror, error, localbnorm, bnorm;
  double tolerance=0.0; //tolerance for convergence. <=0 means do not check

  //main arrays
  double **psi, **zet;
  //temporary versions of main arrays
  double **psitmp, **zettmp;

  //command line arguments
  int scalefactor, numiter;

  double re; // Reynold's number - must be less than 3.7

  //simulation sizes
  int bbase=10;
  int hbase=15;
  int wbase=5;
  int mbase=32;
  int nbase=32;

  int irrotational = 1, checkerr = 0;

  int m,n,lm,b,h,w;
  int iter;
  int i,j;

  /**
   *
   * The code is modified here.
   * These variables are initialised in every rank.
   * tstart, tstop1 and stop2 are the marks for time measurement.
   * ttot_iter, ttot_comp and ttot_comm store the three different types of time we measure
   * **time_iter, **time_comp, **time_comm gather the values from every rank's ttot_iter, ttot_comp and ttot_comm
   *
  **/
  double tstart, tstop1,tstop2,ttot_iter,ttot_comp,ttot_comm,titer;
  double **time_iter,**time_comp,**time_comm;

  //parallelisation parameters
  int rank, size;
  MPI_Comm comm;

  if (tolerance > 0) {checkerr=1;}

  comm=MPI_COMM_WORLD;

  MPI_Init(&argc,&argv);

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  /**
   *
   * The code is modified here.
   * time_iter, time_comp and time_comm are allocated the memory of which size is 5000x432 dynamically
   * 5000 represents the number of iterations
   * 432 is the total cores we used. If users change the number of cores, this value need to be altered.
   *
  **/

  time_iter=(double **)arraymalloc2d(5000,432,sizeof(double));
  time_comp=(double **)arraymalloc2d(5000,432,sizeof(double));
  time_comm=(double **)arraymalloc2d(5000,432,sizeof(double));


  //check command line parameters and parse them

  if (argc <3|| argc >4)
    {
      if (rank == 0) printf("Usage: cfd <scale> <numiter> [reynolds]\n");
      MPI_Finalize();
      return 0;
    }

  if (rank == 0)
    {
      scalefactor=atoi(argv[1]);
      numiter=atoi(argv[2]);

      if (argc == 4)
	{
	  re=atof(argv[3]);
	  irrotational=0;
	}
      else
	{
	  re=-1.0;
	}

  /**
   *
   * We modified the code here.
   * The useless output is removed.
   *
  **/

     /* if(!checkerr)
	{
	  printf("Scale Factor = %i, iterations = %i\n",scalefactor, numiter);
	}
      else
	{
	  printf("Scale Factor = %i, iterations = %i, tolerance= %g\n",scalefactor,numiter,tolerance);
	}

      if (irrotational)
	{
	  printf("Irrotational flow\n");
	}
      else
	{
	  printf("Reynolds number = %f\n",re);
	}*/
    }


  //broadcast runtime params to other processors
  MPI_Bcast(&scalefactor,1,MPI_INT,0,comm);
  MPI_Bcast(&numiter,1,MPI_INT,0,comm);
  MPI_Bcast(&re,1,MPI_DOUBLE,0,comm);
  MPI_Bcast(&irrotational,1,MPI_INT,0,comm);

  //Calculate b, h & w and m & n
  b = bbase*scalefactor;
  h = hbase*scalefactor;
  w = wbase*scalefactor;
  m = mbase*scalefactor;
  n = nbase*scalefactor;

  re = re / (double)scalefactor;

  //calculate local size
  lm = m/size;

  //consistency check
  if (size*lm != m)
    {
      if (rank == 0)
	{
	  printf("ERROR: m=%d does not divide onto %d processes\n", m, size);
	}
      MPI_Finalize();
      return -1;
    }

  /**
   *
   * We modified the code here.
   * The useless output is removed.
   *
  **/
  /*if (rank == 0)
    {
      printf("Running CFD on %d x %d grid using %d process(es) \n",m,n,size);
    }*/

  //allocate arrays

  psi    = (double **) arraymalloc2d(lm+2,n+2,sizeof(double));
  psitmp = (double **) arraymalloc2d(lm+2,n+2,sizeof(double));

  //zero the psi array
  for (i=0;i<lm+2;i++)
    {
      for(j=0;j<n+2;j++)
	{
	  psi[i][j]=0.;
	}
    }

  if (!irrotational)
    {
      //allocate arrays

      zet =   (double **) arraymalloc2d(lm+2,n+2,sizeof(double));
      zettmp =(double **) arraymalloc2d(lm+2,n+2,sizeof(double));

      //zero the zeta array

      for (i=0;i<lm+2;i++)
	{
	  for(j=0;j<n+2;j++)
	    {
	      zet[i][j]=0.;
	    }
	}
    }

  //set the psi boundary conditions

  boundarypsi(psi,lm,n,b,h,w,comm);

  //compute normalisation factor for error

  localbnorm=0.;

  for (i=0;i<lm+2;i++)
    {
      for (j=0;j<n+2;j++)
	{
	  localbnorm += psi[i][j]*psi[i][j];
	}
    }

  //boundary swap of psi

  haloswap(psi,lm,n,comm);

  if (!irrotational)
    {
      //update zeta BCs that depend on psi
      boundaryzet(zet,psi,lm,n,comm);

      //update normalisation

      for (i=0;i<lm+2;i++)
	{
	  for (j=0;j<n+2;j++)
	    {
	      localbnorm += zet[i][j]*zet[i][j];
	    }
	}

      //boundary swap zeta
      haloswap(zet,lm,n,comm);
    }

  //get global bnorm
  MPI_Allreduce(&localbnorm,&bnorm,1,MPI_DOUBLE,MPI_SUM,comm);

  bnorm=sqrt(bnorm);

  //begin iterative Jacobi loop
  /**
   *
   * We modified the code here.
   * The useless output is removed.
   *
  **/
  /*if (rank == 0)
    {
      printf("\nStarting main loop...\n\n");
    }*/
  /**
   *
   * We modified the code here.
   * The location of the timer is changed.
   *
  **/
  //barrier for accurate timing - not needed for correctness

  //MPI_Barrier(comm);

  //tstart=MPI_Wtime();

  for(iter=1;iter<=numiter;iter++)
    {
      /**
       *
       * We modified the code here.
       * This the start of measuring time.
       *
      **/
      //barrier for accurate timing
      MPI_Barrier(comm);
      tstart=MPI_Wtime();

      //calculate psi for next iteration

      if (irrotational)
	{
	  jacobistep(psitmp,psi,lm,n);
	}
      else
	{
	  jacobistepvort(zettmp,psitmp,zet,psi,lm,n,re);
	}

      //calculate current error if required

      if (checkerr || iter == numiter)
	{
	  localerror = deltasq(psitmp,psi,lm,n);

	  if(!irrotational)
	    {
	      localerror += deltasq(zettmp,zet,lm,n);
	    }

	  MPI_Allreduce(&localerror,&error,1,MPI_DOUBLE,MPI_SUM,comm);
	  error=sqrt(error);
	  error=error/bnorm;
	}

      //quit early if we have reached required tolerance

      if (checkerr)
	{
	  if (error < tolerance)
	    {
	      /*if (rank == 0)
		{
		  printf("Converged on iteration %d\n",iter);
		}*/
	      break;
	    }
	}

      //copy back

      for(i=1;i<=lm;i++)
	{
	  for(j=1;j<=n;j++)
	    {
	      psi[i][j]=psitmp[i][j];
	    }
	}

      if (!irrotational)
	{
	  for(i=1;i<=lm;i++)
	    {
	      for(j=1;j<=n;j++)
		{
		  zet[i][j]=zettmp[i][j];
		}
	    }
	}

     /**
       *
       * We modified the code here.
       * An end of the time measurement is added here.
       * The period between this end and the start is the computation time.
       *
      **/
      tstop1 = MPI_Wtime();

      //do a boundary swap

      haloswap(psi,lm,n,comm);

      if (!irrotational)
	{
	  haloswap(zet,lm,n,comm);
	  boundaryzet(zet,psi,lm,n,comm);
	}

      //print loop information
      /**
       *
       * Useless output is removed.
       *
       */
      /*if(iter%printfreq == 0)
	{
	  if (rank==0)
	    {
	      if (!checkerr)
		{
		  printf("Completed iteration %d\n",iter);
		}
	      else
		{
		  printf("Completed iteration %d, error = %g\n",iter,error);
		}
	    }
	}*/
      /**
       *
       * The other end is added here. The time between this end the one before is the communication time.
       * The time before this end and the start is the runtime of iteration.
       * These three types of time are stored in the corresponding variable.
       * The values in the ttot_iter, ttot_comp and ttot_comm from different ranks are gathered in rank 0.
       *
      **/
      tstop2=MPI_Wtime();

      ttot_iter = tstop2 - tstart;
      ttot_comp = tstop1 - tstart;
      ttot_comm = tstop2 - tstop1;
      MPI_Gather(&ttot_iter,1,MPI_DOUBLE,time_iter[iter-1],1,MPI_DOUBLE,0,comm);
      MPI_Gather(&ttot_comp,1,MPI_DOUBLE,time_comp[iter-1],1,MPI_DOUBLE,0,comm);
      MPI_Gather(&ttot_comm,1,MPI_DOUBLE,time_comm[iter-1],1,MPI_DOUBLE,0,comm);
    }

  if (iter > numiter) iter=numiter;

  /**
   *
   * The useless output is removed.
   *
  **/
  //tstop=MPI_Wtime();

  //ttot=tstop-tstart;
  //titer=ttot/(double)iter;

  //MPI_Gather(&ttot,1,MPI_DOUBLE,time,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

  /*
  MPI_Barrier(comm);
  tstop=MPI_Wtime();

  ttot=tstop-tstart;
  titer=ttot/(double)iter;
  */



  /**
   *
   * The data have been created.
   * It outputs the data in some order.
   * The first 5000 rows is the runtime of each iteration.
   * Next 5000 rows is computation time.
   * The last 5000 rows is communication time.
   *
  **/

  if (rank == 0)
    {
      for(i=0;i<5000;i++)
        {
          for(j=0;j<432;j++)
            {
              printf("%f  ",time_iter[i][j]);
            }
	  printf("\n");
         }
      for(i=0;i<5000;i++)
        {
          for(j=0;j<432;j++)
            {
              printf("%f  ",time_comp[i][j]);
            }
          printf("\n");
         }
      for(i=0;i<5000;i++)
        {
          for(j=0;j<432;j++)
            {
              printf("%f  ",time_comm[i][j]);
            }
          printf("\n");
         }
    }

  /**
   *
   * Remove useless data and free the dynamically allocated arrays.
   *
  **/

  //output results

  //writedatafiles(psi,lm,n, scalefactor,comm);

  //if (rank == 0) writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  free(psi);
  free(psitmp);
  free(time_iter);
  free(time_comm);
  free(time_comp);

  if (!irrotational)
    {
      free(zet);
      free(zettmp);
    }

  MPI_Finalize();

  /*if (rank == 0)
    {
      printf("... finished\n");
    }*/

  return 0;
}
