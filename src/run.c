/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
#include "brains.h"

/*! \file run.c
 *  \brief setup and run the program. 
 */

/*!
 * This function setups and runs the program.
 */
void begin_run()
{
  /* Velocity unit 
   * This unit should be determined in the beginning because the data velocity
   * needs to be converted using this unit in read.c. 
   */
  VelUnit = sqrt( GRAVITY * 1.0e6 * SOLAR_MASS / CM_PER_LD ) / 1.0e5; 

#ifdef SA
  PhaseFactor = (1.0/360.0 * CM_PER_PC/CM_PER_LD);
#endif

  /* dimensionless speed of light */
  C_Unit = C/1.0e5/VelUnit;

  /* read parameter file */
  read_parset();
  if(thistask == roottask)
  {
    fprint_param();
  }

  if(parset.flag_dim != -2) /* if not randomly create mock data */
  {
    /* read data files */
    read_data();
#ifndef SA
    /* scale continuum and line to an order of unity */
    scale_con_line();
#else
    if(parset.flag_dim != 3)
    {
      scale_con_line();
    }
#endif
  }

  /* initialization */
  init();
  
  /* now run dnest and reconstruct the model. */
  MPI_Barrier(MPI_COMM_WORLD);


  if(parset.flag_dim == -1 || parset.flag_dim == -2) /* only simulation */
  {
    /* first determine the best variability parameters in continuum */
    if(parset.flag_dim == -1)
      reconstruct_con();
    
    sim();
  }

  if(parset.flag_dim == 0) /* only continuum */
  {
    reconstruct_con();
  }

  if(parset.flag_dim == 1) /* 1d line */
  {
    if(parset.flag_postprc == 0)
    {
      reconstruct_con();
    }
  	reconstruct_line1d();
  }

  if(parset.flag_dim == 2) /* 2d line */
  {
    if(parset.flag_postprc == 0)
    {
      reconstruct_con();
    }
    reconstruct_line2d();
  }

#ifdef SA
  if(parset.flag_dim == 3) /* SA */
  {
    reconstruct_sa();
  }

  if(parset.flag_dim == 4) /* SA + RM */
  {
    if(parset.flag_postprc == 0)
    {
      reconstruct_con();
    }
    reconstruct_sa1d();
  }

  if(parset.flag_dim == 5) /* SA + RM */
  {
    if(parset.flag_postprc == 0)
    {
      reconstruct_con();
    }
    reconstruct_sa2d();
  }
#endif

  return;
}

/*!
 * This function frees the memory and ends the run.
 */
void end_run()
{
  free_memory_data();
  free_memory();
}