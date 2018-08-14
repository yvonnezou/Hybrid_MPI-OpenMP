#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <mpi.h>
#include <omp.h>

#include "arraymalloc.h"
#include "boundary.h"
#include "jacobi.h"
#include "cfdio.h"

int main(int argc, char **argv)
{
  int printfreq=1000; //output frequency
  double localerror, error, localbnorm,bnorm;
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

  int required,provided;
  required = MPI_THREAD_MULTIPLE;

  double tstart, tstop, ttot, titer;

  double **time_iter, **time_comm, **time_comp, **time_iter_thread, **time_comm_thread, **time_comp_thread;

  //initial time array
  time_iter    = (double **) arraymalloc2d(5000,384,sizeof(double));
  time_comm    = (double **) arraymalloc2d(5000,384,sizeof(double));
  time_comp    = (double **) arraymalloc2d(5000,384,sizeof(double));
  time_iter_thread = (double **)arraymalloc2d(5000,12,sizeof(double));
  time_comm_thread = (double **)arraymalloc2d(5000,12,sizeof(double));
  time_comp_thread = (double **)arraymalloc2d(5000,12,sizeof(double));

/*
  for(i=0;i<5000;i++){
     for(j=0;j<384;j++){
        time_iter[i][j] = 0;
	time_comm[i][j] = 0;
        time_comp[i][j] = 0;
     }
  }
*/

  //parallelisation parameters
  int rank, size;
  //int nthreads;			//number of threads & thread id
  MPI_Comm comm;


  //do we stop because of tolerance?
  if (tolerance > 0) {checkerr=1;}

  comm=MPI_COMM_WORLD;

  MPI_Init_thread(&argc,&argv,required,&provided);

  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

/*
#pragma omp parallel shared(nthreads)
  nthreads = omp_get_num_threads();
  //ithread = OMP_GET_THREAD_NUM();
*/
  
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
	
      /*if(!checkerr)
	{
	  printf("Scale Factor = %i, iterations = %i\n",scalefactor, numiter);
	}
      else
	{
	  printf("Scale Factor = %i, iterations = %i, tolerance= %g\n",scalefactor,numiter,tolerance);
	}
      */

      /*if (irrotational)
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

  /*if (rank == 0)
    {
      printf("Running CFD on %d x %d grid using %d process(es)\n",m,n,size);
    }
  */

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

  localbnorm=0.0;

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

  //if (rank == 0)
  //  {
  //    printf("\nStarting main loop...\n\n");
  //  }

  MPI_Barrier(comm);

  //tstart=MPI_Wtime();
  
  //OpenMP region
#pragma omp parallel default (none) shared(time_iter_thread,time_comm_thread,time_comp_thread,zettmp,psitmp,zet,psi,n,numiter,irrotational,lm,re,comm,printfreq,rank,checkerr) private(i,j)
  {
  int nthreads = omp_get_num_threads();
  int myid = omp_get_thread_num();
  double tstart_thread,tstop1_thread,tstop2_thread;
  int iter;

  //The size of each chunk in every thread
  int ln = n/nthreads;
  
  if (n%nthreads != 0)
    {
      printf("unequal/n");
    } 

  for(iter=0;iter<numiter;iter++)
    {     

      tstart_thread = omp_get_wtime();

      if (irrotational)
	{
	  jacobistep(psitmp,psi,lm,ln,myid);
	}
      else
	{
	  jacobistepvort(zettmp,psitmp,zet,psi,lm,ln,re,myid,ln);
	}
      //copy back ***********
      for(i=1;i<=lm;i++)
	{
	  for(j=1+myid*ln;j<=(1+myid)*ln;j++)
	    {
	      psi[i][j]=psitmp[i][j];
	    }
	}

      if (!irrotational)
	{
	  for(i=1;i<=lm;i++)
	    {
	      for(j=1+myid*ln;j<1+(1+myid)*ln;j++)
		{
		  zet[i][j]=zettmp[i][j];
		}
	    }
	}

      tstop1_thread = omp_get_wtime();
      
      //do a boundary swap ************
      haloswap_thread(psi,lm,ln,comm,myid);

      tstop2_thread = omp_get_wtime();
      
      //timer in every thread

      time_iter_thread[iter][myid] = tstop2_thread - tstart_thread;
      time_comp_thread[iter][myid] = tstop1_thread - tstart_thread;
      time_comm_thread[iter][myid] = tstop2_thread - tstop1_thread;
      }

      //print loop information

      /*if(iter%printfreq == 0)
	{
	  printf("Rank %d Thread %d Completed iteration %d\n",rank,myid,iter);
	}
      */
    

  }

  //if (iter > numiter) iter=numiter;

  MPI_Barrier(comm);
  //tstop=MPI_Wtime();

  //ttot=tstop-tstart;
  //titer=ttot/(double)numiter;


  //print out some stats
  /*if (rank == 0)
    {
      printf("\n... finished\n");
      //printf("After %d iterations, the error is %g\n",iter,error);
      printf("Time for %d iterations was %g seconds\n",numiter,ttot);
      printf("Each iteration took %g seconds\n",titer);
    }
  */

  //print time in the thread
  //printf("rank = %d\n",rank);
  for(i=0;i<5000;i++) {
     MPI_Gather(time_iter_thread[i],12,MPI_DOUBLE,time_iter[i],12,MPI_DOUBLE,0,comm);
     MPI_Gather(time_comp_thread[i],12,MPI_DOUBLE,time_comp[i],12,MPI_DOUBLE,0,comm);
     MPI_Gather(time_comm_thread[i],12,MPI_DOUBLE,time_comm[i],12,MPI_DOUBLE,0,comm);
  }

  //print runtime
  if(rank == 0) {
     for(i=0;i<5000;i++) {
        for(j=0;j<384;j++) {
	    printf("%f  ",time_iter[i][j]);
        }
        printf("\n");
     }

     for(i=0;i<5000;i++) {
        for(j=0;j<384;j++) {
              printf("%f  ",time_comp[i][j]);
         }
         printf("\n");
     }

     for(i=0;i<5000;i++) {
        for(j=0;j<384;j++) {
              printf("%f  ",time_comm[i][j]);
         }
         printf("\n");
     }
   }

  //output results

  //writedatafiles(psi,lm,n, scalefactor,comm);

  //if (rank == 0) writeplotfile(m,n,scalefactor);

  //free un-needed arrays
  free(psi);
  free(psitmp);
  free(time_iter);
  free(time_comp);
  free(time_comm);
  free(time_comm_thread);
  free(time_comp_thread);
  free(time_iter_thread);

  if (!irrotational)
    {
      free(zet);
      free(zettmp);
    }

  MPI_Finalize();
  return 0;
}
