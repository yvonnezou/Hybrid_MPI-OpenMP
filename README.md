# Characterising performance noise on HPC systems
---
This project consists of a benchmark in three different style, pure-MPI programme, hybrid programme in single style and hybrid programme in multiple style. The benchmark is a paralleled Fluid Dynamics programme using C language, MPI and OpenMP. These programmes are designed to measure the noise on the HPC system ARCHER and Cirrus. It also includes the data processing scripts and the results from these two HPC systems.

### Table of contents
---
1. File description
2. Usage

### File description
---
The three different kinds of programme are stored in the folder of `Pure_MPI programme`, `Hybrid programme in single style` and `Hybrid programme in multiple style` respectively. The programme in single style does not suit the experiments. The other two programmes have been put more emphasis on.

In each folder, there are two folders named `Archer` and `Cirrus` storing corresponding programmes. Inside, there are several files:
 - arraymalloc.c: Allocate the arrays dynamically
 - boundary.c: Swap the halos
 - cfd.c: Main file
 - cfdio.c: Output files
 - jacobi.c: Do the Jacobi Algorithm calculations
 - cfd.pbs/t2cfd.pbs: The script to submit the programme to backend
 - runsplit.sh: Process the data
 - run.sh: Submit the programme to backend 100 times

In the Data analytics and results folder, the folders, `3D images`, `combined images` and `detailed images` store different kind of output. In `combined image`, the maximum, minimum and average value are in the same figure. In `detailed image`, they are separate. There also two Python script inside:
 - 3d.py: Plot 3D images
 - plot.py: Plot bar charts

### Usage
---
1. Download the package and upload to ARCHER or Cirrus by `scp` command.
2. Enter the corresponding folder and `make` the programme. On Cirrus, it is necessary to load module intel-compilers-17 and intel-mpi-17 first.
    ```
    $ module load intel-compilers-17
    $ module load intel-mpi-17
    ```
3. Submit it to backend. 
    On Archer:
    ```
    $ qsub cfd.pbs
    ```
    On Cirrus:
    ```
    $ qsub t2cfd.pbs/t2cfd_hpempt.pbs/t2cfd_intel.pbs
    ````
4. Process the data
    ```
    $ sh runsplit.sh
    ```
5. Download the txt data file to local by using `scp`
6. Store the data and Python script in the same path and tun the Python script to get results.
