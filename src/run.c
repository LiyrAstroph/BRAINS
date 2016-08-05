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
  read_parset();
  read_data();
  scale_con_line();

  init();

  if(parset.flag_only_recon)
  {
    reconstruct_con();
  }
}

void end_run()
{
  free_memory_data();
}