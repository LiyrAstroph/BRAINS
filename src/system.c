/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file system.c
 *  \brief get time used by the program.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <signal.h>
#include <mpi.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

/*!
 * This function provides the clock time of the system.
 */
double second()
{
#ifdef WALLCLOCK
	return MPI_Wtime();
#else
	return ((double) clock()) /CLOCKS_PER_SEC;
#endif
}

/*!
 *  This function calculates the time difference.
 */
double timediff(double t0, double t1)
{
	double dt;
	dt = t1 -t0;
	if(dt<0)
	{
#ifdef WALLCLOCK
	dt =0;
#else   
	dt = t1 + pow(2, sizeof(clock_t)*8) / CLOCKS_PER_SEC - t0;  /* 1 byte = 8 bits */
#endif
  }
	return dt;
}

/*!
 *  This function converts clock time into hms.
 */
void get_hms(double dt, int *h, int *m, double *s)
{
	*h = (int) floor(dt/3600);
	*m = (int) floor((dt - (*h)*3600)/60);
  *s = dt - (*h)*3600 - (*m)*60;
  return;
}