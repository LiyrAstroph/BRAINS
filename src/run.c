/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "allvars.h"
#include "proto.h"

void begin_run()
{
  /* Velocity unit */
  VelUnit = sqrt( GRAVITY * 1.0e6 * SOLAR_MASS / CM_PER_LD ) / 1.0e5;

  read_parset();
  read_data();
  scale_con_line();

  init();
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(parset.flag_dim == 0)
  {
    reconstruct_con();
  }

  if(parset.flag_dim == 1)
  {
  	reconstruct_line1d();
  }

  if(parset.flag_dim == 2)
  {
    reconstruct_line2d();
  }
}

void end_run()
{
  free_memory_data();
  free_memory();
}