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

#define NONE                 "\e[0m"
#define BLACK                "\e[0;30m"
#define L_BLACK              "\e[1;30m"
#define RED                  "\e[0;31m"
#define L_RED                "\e[1;31m"
#define GREEN                "\e[0;32m"
#define L_GREEN              "\e[1;32m"
#define BROWN                "\e[0;33m"
#define YELLOW               "\e[1;33m"
#define BLUE                 "\e[0;34m"
#define L_BLUE               "\e[1;34m"
#define PURPLE               "\e[0;35m"
#define L_PURPLE             "\e[1;35m"
#define CYAN                 "\e[0;36m"
#define L_CYAN               "\e[1;36m"
#define GRAY                 "\e[0;37m"
#define WHITE                "\e[1;37m"

#define BOLD                 "\e[1m"
#define UNDERLINE            "\e[4m"
#define BLINK                "\e[5m"
#define REVERSE              "\e[7m"
#define HIDE                 "\e[8m"
#define CLEAR                "\e[2J"
#define CLRLINE              "\r\e[K" //or "\e[1K\r"

void print_help()
{
  printf("\n\n");
  printf(BOLD "NAMES\n" NONE);
  printf("\tBRAINS --- an MPI-based code for analyzing reverberation mapping data.\n" NONE);
  printf("\n");
  
  printf(BOLD "SYNOPSIS\n" NONE);
  printf("\t" BOLD "brains [FILE] [OPTION]\n" NONE);
  printf("\t" BOLD "mpiexec [MPI_OPTION] brains [FILE] [OPTION]\n" NONE);
  printf("\n");

  printf(BOLD "OPTIONS\n" NONE);
  printf("\tMPI_OPTIONS\n");
  printf("\t\toptions passed to MPI.\n");

  printf("\tFILE\n");
  printf("\t\tparameter file.\n");

  printf("\t" PURPLE BOLD "-h or --help\n" NONE);
  printf("\t\tprint help information.\n");

  printf("\t" PURPLE BOLD "-p or --post_proc\n" NONE);
  printf("\t\tonly do posterior processing.\n");

  printf("\t" PURPLE BOLD "-r or --restart\n" NONE);
  printf("\t\trestart from the backup.\n");

  printf("\t" PURPLE BOLD "-t or --temperature or --temp\n" NONE);
  printf("\t\tspecify tempering temperature in posterior processing.\n");

  printf("\t" PURPLE BOLD "-s or --seed\n" NONE);
  printf("\t\tset a seed for the random number generator.\n");

  printf("\t" PURPLE BOLD "-c or --recalc_info\n" NONE);
  printf("\t\tonly do posterior processing, but recalculate the posterior sample information.\n");

  printf("\t" PURPLE BOLD "-e or --exam_prior\n" NONE);
  printf("\t\texamine the priors.\n");

  printf("\t" PURPLE BOLD "-v or --version\n" NONE);
  printf("\t\tprint version.\n");

  printf("\t" PURPLE BOLD "-n or --para_name\n" NONE);
  printf("\t\tprint parameter name.\n");

  printf("\n\n");
  return;
}