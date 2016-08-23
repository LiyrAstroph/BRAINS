/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>
 
#include "dnest_line2d.h"
#include "allvars.h"
#include "proto.h"

void reconstruct_line2d()
{
  reconstruct_line2d_init();
  
  dnest_line2d(0, " ");

  if(thistask == roottask)
  {
    calculate_con_from_model(best_model_line2d + 11 *sizeof(double));
    FILE *fp;
    char fname[200];
    int i, j;

    sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tcon[i], Fcon[i] / con_scale);
    }
    fclose(fp);

    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

    // recovered line2d at data points
    which_parameter_update = 1;
    Trans2D_at_veldata = Trans2D_at_veldata_particles[0];
    transfun_2d_cloud_direct(best_model_line2d, Vline_data, Trans2D_at_veldata, 
    	                                        n_vel_data, parset.flag_save_clouds);
    calculate_line2d_from_blrmodel(best_model_line2d, Tline_data, Vline_data, Trans2D_at_veldata,
    	                                            Fline2d_at_data, n_line_data, n_vel_data);

    sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_data_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<n_line_data; i++)
    {
      for(j=0; j<n_vel_data; j++)
      {
      	fprintf(fp, "%f %f %f\n", Vline_data[j]*VelUnit, Tline_data[i],  Fline2d_at_data[i*n_vel_data + j] / line_scale);
      }

      fprintf(fp, "\n");
    }
    fclose(fp);

    // recovered line2d at specified points
    transfun_2d_cloud_direct(best_model_line2d, TransV, Trans2D, parset.n_vel_recon, 0);
    calculate_line2d_from_blrmodel(best_model_line2d, Tline, TransV, 
    	    Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

    sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_line_recon; i++)
    {
      for(j=0; j<parset.n_vel_recon; j++)
      {
      	fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, Tline[i],  Fline2d[i*parset.n_vel_recon + j] / line_scale);
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

    for(i=0; i<parset.n_line_recon; i++)
    {
      for(j=0; j<parset.n_vel_recon; j++)
      {
        fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, TransTau[i], Trans2D[i*parset.n_vel_recon + j]);
      }

      fprintf(fp, "\n");
    }
    fclose(fp);
    
  }

  reconstruct_line2d_end();
}

void reconstruct_line2d_init()
{
  TransTau = malloc(parset.n_tau * sizeof(double));
  //Trans2D_at_veldata = malloc(parset.n_tau * n_vel_data * sizeof(double));


  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans2D = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));
  
  Fline2d_at_data = malloc(n_line_data * n_vel_data * sizeof(double));

  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  int i;
  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 20);
  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 20);
  double dT = (Tline_max - Tline_min)/(n_line_data - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }

  dTransTau = (parset.tau_max_set - parset.tau_min_set)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = parset.tau_min_set + dTransTau * i;
  }

  double vel_max_set = Vline_data[n_vel_data -1], vel_min_set = Vline_data[0];
  double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon);

  for(i=0; i<parset.n_vel_recon; i++)
  {
  	TransV[i] = vel_min_set + dVel*i;
  }

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS2D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  Trans2D_at_veldata_particles = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans2D_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_data * sizeof(double));
  }

  // only record gamma-distribution random number of clouds
  clouds_particles = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    clouds_particles[i] = malloc(parset.n_cloud_per_task * sizeof(double));
  }

  beta_old_particles = malloc(parset.num_particles * sizeof(double));
  for(i=0; i<parset.num_particles; i++)
    beta_old_particles[i] = -1.0;

  return;
}

void reconstruct_line2d_end()
{
  free(Tline);
  //free(Fline2d_at_data);
  free(Fline2d);

  free(TransTau);
  free(TransV);
  free(Trans2D_at_veldata);
  free(Trans2D);

  int i;
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans2D_at_veldata_particles[i]);
    free(clouds_particles[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(clouds_particles);
  free(beta_old_particles);
  return;
}

double prob_line2d(void *model)
{
  double prob = 0.0, fcon, var2, dy;
  int i;
  
  calculate_con_from_model(model + 11*sizeof(double));
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

  for(i=0; i<n_con_data; i++)
  {
    fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
    prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - ( 0.5*log(2.0*PI) + log(Fcerrs_data[i]) );
  }
  
  // only update transfer function when BLR model is changed.
  if(which_parameter_update < 11)
  {
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
    transfun_2d_cloud_direct(model, Vline_data, Trans2D_at_veldata, n_vel_data, 0);
    //memcpy(Trans2D_at_veldata_particles[which_particle_update], Trans2D_at_veldata, 
    //  n_vel_data*parset.n_tau * sizeof(double));
  }
  else
  {
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
    //memcpy(Trans2D_at_veldata, Trans2D_at_veldata_particles[which_particle_update], 
    //  n_vel_data*parset.n_tau * sizeof(double));
  }
  calculate_line2d_from_blrmodel(model, Tline_data, Vline_data, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data);

  for(i=0; i<n_line_data*n_vel_data; i++)
  {
    dy = Fline2d_data[i] - Fline2d_at_data[i] ;
    var2 = Flerrs2d_data[i]*Flerrs2d_data[i];
    var2 = 1.0*var2;
    prob += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
  }
  if(isnan(prob))
    prob = -DBL_MAX;
  return prob;
}

