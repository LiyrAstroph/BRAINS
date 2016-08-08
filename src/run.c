/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>

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

  if(parset.flag_only_recon)
  {
    reconstruct_con();
  }

  if(!parset.flag_only_recon && parset.flag_dim == 1)
  {
  	reconstruct_line1d();
  }

  if(!parset.flag_only_recon && parset.flag_dim == 2)
  {
    reconstruct_line2d();
  }
}

void end_run()
{
  free_memory_data();
}