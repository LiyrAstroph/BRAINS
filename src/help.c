/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h> 
#include <unistd.h>
#include <mpi.h>
#include <string.h>

#include "allvars.h"
#include "proto.h"

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

  printf("\t-h\n");
  printf("\t\tprint help information.\n");

  printf("\t-p\n");
  printf("\t\tonly do posterior processing.\n");

  printf("\t-r\n");
  printf("\t\trestart from the backup.\n");

  printf("\t-t\n");
  printf("\t\tspecify tempering temperature in posterior processing.\n");

  printf("\t-c\n");
  printf("\t\tonly do posterior processing, but recalculate the posterior sample information.\n");

  printf("\n\n");
  return;
}