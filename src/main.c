/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

/*! \file main.c
 *  \brief start of the program
 */
 
/*!
 *  This function initializes the MPI communication packages.
 */ 
int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &thistask);
  MPI_Comm_size(MPI_COMM_WORLD, &totaltask);
  MPI_Get_processor_name(proc_name, &namelen);
  
  if(thistask == roottask)
  {
    printf("===============BRAINS==================\n");
    printf("Starts to run...\n");
    printf("%d cores used.\n", totaltask);
  }

  if(argc<2)
  {
    if(thistask == roottask)
    {
      fprintf(stderr, "# Error: No parameter file specified!\n");
      fprintf(stdout, "Ends incorrectly.\n");
    }
    MPI_Finalize();
    return 0;
  }

  strcpy(parset.param_file, argv[1]); /* copy input parameter file */


  begin_run();    /* implementation run */

  end_run();      /* end run */
  
  MPI_Finalize();   /* clean up and finalize MPI */
  if(thistask == roottask)
  {
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
   <b>A</b>nalysis <b>I</b>ntegrated with <b>N</b>ested <b>S</b>ampling.

\section install Compilation 

BRAINS needs the following non-standard libraries for compilation:

- \b MPI - the Message Passing Interface (version 1.0 or higher).

\section howtorun Running the code
 */

