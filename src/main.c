/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <mpi.h>
#include <string.h>

#include "brains.h"

/*! \file main.c
 *  \brief start of the program
 */
 
/*!
 *  This function initializes the MPI communication packages.
 */ 
int main(int argc, char **argv)
{
  double t0=0.0, t1, dt;
  
  /* initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  
  if(thistask == roottask)
  {
    t0 = second();
    printf("===============BRAINS==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }
  
  if(command_line_options(argc, argv) != EXIT_SUCCESS)
  {
    MPI_Finalize();
    return 0;
  }

  if(parset.flag_help == 0)
  {
    begin_run();    /* implementation run */

    end_run();      /* end run */
  }
  
  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
    int ht, mt;
    double st;
    t1 = second();
    dt = timediff(t0, t1);
    get_hms(dt, &ht, &mt, &st);
    printf("Time used: %dh %dm %fs.\n", ht, mt, st);
    printf("Ends successfully.\n");
    printf("===============BRAINS==================\n");
  }
  return 0;
}

/* ----------------------------------------------------------------------
   The rest of this file contains documentation for compiling and
   running the code, in a format appropriate for 'doxygen'.
   ----------------------------------------------------------------------
 */

/*! \mainpage Reference documentation of BRAINS

\author Yan-Rong Li\n
        liyanrong@ihep.ac.cn\n
        Institute of High Energy Physics

\section introduction Introduction

\b BRAINS is an abbreviation for  <b>B</b>LR <b>R</b>everberation-mapping 
   <b>A</b>nalysis <b>I</b>n AGNs with <b>N</b>ested <b>S</b>ampling.

\section install Compilation 

BRAINS needs the following non-standard libraries for compilation:

- \b MPI - the Message Passing Interface (version 1.0 or higher).\n
  The widely used MPI implementation MPICH is recommended, which can 
  be downloaded at http://www-unix.mcs.anl.gov/mpi/mpich/. \n Note that 
  many vendor supplied versions exist.

- \b GSL - the GNU Scientific Library (version 2.2.1 or higher).\n
  This open-source package can be downloaded at http://www.gnu.org/software/gsl.\n
  BRAINS uses this library to generate random numbers, perform FFT
  tranformation, and access to CBLAS library. 

- \b LAPACKE - the C-interface of LAPACK (version 3.6.1 or higher).\n
  This open-source package can be downloaded at http://www.netlib.org/lapack/.\n
  BRAINS uses this library to perform linear algebra computation.
  Note that to LAPACKE is contained in LAPACK source. One needs to switch on
  LAPACKE when compiling LAPACK.

- \b FFTW - a fast Fourier transform library (version 3.0 or higher).\n
  This open-source package can be downloaded at http://www.fftw.org/.\n

- \b CDNest - a package for diffusive nested sampling.\n
  This package can be downloaded at https://github.com/LiyrAstroph/CDNest.\n
  

Note that in Linux system, there are package managers that can install the above 
libraries convienently (except for CDNest). If so, use them. In this case, the libraries usually are 
installed in standard environment path. Otherwise, if any of the above libraries is 
not installed in standard locations on your system, the \ref BRAINS-Makefile 
"Makefile" provided with the code may need slight adjustments.

\section howtorun Running the code

To run the code in a terminal, type:\n
\code{.sh} 
mpiexec -n np ./brains param
\endcode
where np is the number of cores, param is the parameter file.
In order to start the code, a \ref parameterfile "parameterfile" 
and a \ref optionsfile "optionfile" needs to be specified.
 */

/*! \page BRAINS-Makefile  Makefile of BRAINS
 Makefile
 --------

 <table>
 <tr><th> Command                     <th> Note 
 <tr><td>SHELL=/bin/bash              <td> Bash shell
 <tr><td>CC       = mpicc             <td> MPI C Compiler
 <tr><td>OPTIMIZE = -O2 -Wall         <td> Compile option
 <tr><td>#OPTIMIZE += -DDebug         <td> Debug with fixed random seed
 <tr><td>#OPTIMIZE ++ -DSA            <td> Include SA analysis
 <tr><td>
     #SYSTEM="Darwin"\n
      SYSTEM="Linux" \n
     #SYSTEM="Cluster" \n                 <td> Target system
 </table>

 */