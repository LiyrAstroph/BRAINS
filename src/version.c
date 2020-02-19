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
  printf("\n");
  printf("\e[1;35m" "%-14s: %d.%d.%d\n" "\e[0m", "BRAINS Version", BRAINS_MAJOR_VERSION, BRAINS_MINOR_VERSION, BRAINS_PATCH_VERSION);
#ifdef GITVERSION
  printf("\e[1;35m" "%-14s: %s\n" "\e[0m", "git log", GITVERSION);
#endif
#ifdef GITDATE
  printf("\e[1;35m" "%-14s: %s\n" "\e[0m", "git date", GITDATE);
#endif
  printf("\e[1;35m" "%-14s: %s %s\n" "\e[0m", "compiling date", __DATE__, __TIME__);
  printf("\n");
  printf("Yan-Rong Li, liyanrong@ihep.ac.cn\n");
  return;
}