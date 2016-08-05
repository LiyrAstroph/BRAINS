/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void init()
{
  allocate_memory();

  int i;
  double Tcon_min, Tcon_max, dT;
  
  /* set time array for continuum */
  Tcon_min = Tcon_data[0] - fmax(0.1*(Tcon_data[n_con_data -1] - Tcon_data[0]), 20);
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.1*(Tcon_data[n_con_data -1] - Tcon_data[0]), 20);
  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  /* initialize GSL */
  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  gsl_rng_set(gsl_r, time(NULL)+thistask+1350); 

  gsl_acc = gsl_interp_accel_alloc();
  gsl_linear = gsl_interp_alloc(gsl_interp_linear, parset.n_con_recon);
}

void allocate_memory()
{
  Tcon = malloc(parset.n_con_recon * sizeof(double));
  Fcon = malloc(parset.n_con_recon * sizeof(double));
  Fcerrs = malloc(parset.n_con_recon * sizeof(double));

  PSmat = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));
}

void free_memory()
{
  free(Tcon);
  free(Fcon);
  free(Fcerrs);

  free(PSmat);
}

void scale_con_line()
{
  int i, j;
  double ave_con, ave_line;
  
  con_scale = 1.0;
  line_scale = 1.0;

  ave_con = 0.0;
  ave_line = 0.0;
  for(i=0; i<n_con_data; i++)
  {
    ave_con += Fcon_data[i];
  }

  ave_con /= n_con_data;
  con_scale = 1.0/ave_con;

  for(i=0; i<n_con_data; i++)
  {
    Fcon_data[i] *=con_scale;
    Fcerrs_data[i] *=con_scale;
  }

  printf("con scale: %e\t%e\n", con_scale, ave_con);
  
  if(parset.flag_only_recon)
  {
    return;
  }

  for(i=0; i<n_line_data; i++)
  {
    ave_line += Fline_data[i];
  }
  ave_line /=n_line_data;

  line_scale = 1.0/ave_line;
  
  printf("line scale: %e\t%e\n", line_scale, ave_line);

  for(i=0; i<n_line_data; i++)
  {
    Fline_data[i] *= line_scale;
    Flerrs_data[i] *= line_scale;

    if(parset.flag_dim==2)
      for(j=0; j<n_vel_data; j++)
      {
        Fline2d_data[i][j] *= line_scale;
        Flerrs2d_data[i][j] *= line_scale;
      } 
  }
}