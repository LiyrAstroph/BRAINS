/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>

#include "dnest_line1d.h"
#include "allvars.h"
#include "proto.h"

void reconstruct_line1d()
{
  reconstruct_line1d_init();

  dnest_line1d(0, " ");
}

void reconstruct_line1d_init()
{

}