/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>

#include "dnest_line2d.h"
#include "allvars.h"
#include "proto.h"

void reconstruct_line2d()
{
  reconstruct_line2d_init();

  dnest_line2d(0, " ");
}

void reconstruct_line2d_init()
{

}