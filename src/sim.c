/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */


/*!
 * \file sim.c
 * \brief generate mocked 2d data.
 */
#include <stdio.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

void *model;

void sim()
{
  if(thistask != roottask)
    return;

  FILE *fp;
  char fname[200];
  int i, j;

  printf("%f\n", line_scale);

  sim_init();
  
  double *pm = (double *)model;
  //model1
  /*pm[0] = log(1.0);
  pm[1] = 0.0;
  pm[2] = log(4.0);
  pm[3] = 0.9;
  pm[4] = 0.2;
  pm[5] = 20.0;
  pm[6] = 40.0;
  pm[7] = 0.0;
  pm[8] = log(3.0);
  pm[9] = 0.1;
  pm[10] = 0.5;*/

  //model5
  pm[0] = log(1.0);
  pm[1] = 0.0;
  pm[2] = log(4.0);
  pm[3] = 0.5;
  pm[4] = 2.0;
  pm[5] = 1.5;
  pm[6] = 20.0;
  pm[7] = 40.0;
  pm[8] = 0.0;
  pm[9] = 1.0;
  pm[10] = 1.0;
  pm[11] = log(2.0);
  pm[12] = 0.5;
  pm[13] = 0.5;
  pm[14] = 0.01;
  pm[15] = 0.1;
  pm[16] = 0.01;
  pm[17] = 0.1;
  pm[18] = 0.0;
  pm[19] = log(1.0e-10);


  smooth_init(parset.n_vel_recon, TransV);
  
  reconstruct_con_from_varmodel(0.26, 74.0, 1.0, 0.0);
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_con_full.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    fprintf(fp, "%f %f\n", Tcon[i], Fcon[i]/con_scale);
  }
  fclose(fp);
  
  transfun_1d_cloud_direct(model, 0);
  calculate_line_from_blrmodel(model, Tline, Fline, parset.n_line_recon);

  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_hb.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  for(i=0; i<parset.n_line_recon; i++)
  {
    fprintf(fp, "%f %f %f\n", Tline[i], Fline[i] + + gsl_ran_gaussian(gsl_r, 0.01), 0.01);
  }

  transfun_2d_cloud_direct(model, TransV, Trans2D, parset.n_vel_recon, 0);
  calculate_line2d_from_blrmodel(model, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_hb2d.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  fprintf(fp, "# %d %d\n", parset.n_line_recon, parset.n_vel_recon);
  for(i=0; i<parset.n_line_recon; i++)
  {
    for(j=0; j<parset.n_vel_recon; j++)
    {
      fprintf(fp, "%f %f %f %f\n", TransV[j]*VelUnit, Tline[i],  
        (Fline2d[i*parset.n_vel_recon + j] + gsl_ran_gaussian(gsl_r, 0.01)), 0.01);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);
  
  // output 2d transfer function
  sprintf(fname, "%s/%s", parset.file_dir, parset.tran2d_out_file);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  fprintf(fp, "# %d %d\n", parset.n_tau, parset.n_vel_recon);
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<parset.n_vel_recon; j++)
    {
      fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, TransTau[i], Trans2D[i*parset.n_vel_recon + j]);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);

  smooth_end();
  sim_end();
}

void sim_init()
{
  int i;
  double dT, Tspan;

  switch(parset.flag_blrmodel)
  {
    case 1:
      num_params_blr = 12;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model1;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model1;
      break;
    case 2:
      num_params_blr = 12;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model1;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model2;
      break;
    case 3:
      num_params_blr = 12;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model3;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model3;
      break;
    case 4:
      num_params_blr = 12;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model3;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model4;
      break;
    case 5:
      num_params_blr = 20;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model5;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model5;
      break;
    default:
      num_params_blr = 12;
      transfun_1d_cloud_direct = transfun_1d_cloud_direct_model1;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model1;
      break;
  }

  parset.num_particles = 1;
  which_particle_update = 0;
  force_update = 1;
  which_parameter_update = -1;
  
  num_params_var = 4 + parset.flag_trend;
  num_params = num_params_blr + num_params_var + parset.n_con_recon;

  model = malloc(num_params * sizeof(double));

  Fcon = malloc(parset.n_con_recon * sizeof(double));

  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  
  /* set time array for continuum */
  Tcon_min = Tcon_data[0] - fmax(0.05*Tspan, parset.tau_max_set);
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 10.0);
  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  TransTau = malloc(parset.n_tau * sizeof(double));
  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans1D = malloc(parset.n_tau * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  Tline_min = Tcon_min + fmax(0.05*Tspan, parset.tau_max_set);
  Tline_max = Tcon_max;

  dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }
  
  dTransTau = (parset.tau_max_set - parset.tau_min_set)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = parset.tau_min_set + dTransTau * i;
  }

  double vel_max_set, vel_min_set;
  vel_max_set = 3000.0/VelUnit;
  vel_min_set = - vel_max_set;
  double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

  for(i=0; i<parset.n_vel_recon; i++)
  {
    TransV[i] = vel_min_set + dVel*i;
  }
  
  clouds_particles = malloc(parset.num_particles * sizeof(double *));
  clouds_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    clouds_particles[i] = malloc(parset.n_cloud_per_task * sizeof(double));
    clouds_particles_perturb[i] = malloc(parset.n_cloud_per_task * sizeof(double));
  }

  con_q = malloc(nq * sizeof(double));
  return;
}

void sim_end()
{
  int i;
  free(model);
  free(Fcon);

  free(TransTau);
  free(TransV);
  free(Tline);
  free(Fline);
  free(Fline2d);
  free(Trans2D);
  free(Trans1D);
  for(i=0; i<parset.num_particles; i++)
  {
    free(clouds_particles[i]);
    free(clouds_particles_perturb[i]);
  }
  free(clouds_particles);
  free(clouds_particles_perturb);

  free(con_q);
}
