/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>

#include "version.h"
#include "proto.h"

void print_version()
{
  printf("\e[1;35m" "BRAINS Version: %d.%d.%d\n" "\e[0m", BRAINS_MAJOR_VERSION, BRAINS_MINOR_VERSION, BRAINS_PATCH_VERSION);
  printf("Yan-Rong Li, liyanrong@ihep.ac.cn\n");
  return;
}