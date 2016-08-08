/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <mpi.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"
 
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

  strcpy(parset.param_file, argv[1]);

  begin_run();

  end_run();
  
  MPI_Finalize();
  if(thistask == roottask)
  {
    printf("Ends successfully.\n");
    printf("===============BRAINS==================\n");
  }
  return 0;
}