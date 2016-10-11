/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>

#include "dnestvars.h"

#include "dnest_line2d.h"
#include "allvars.h"
#include "proto.h"

void *best_model_line2d, *best_model_std_line2d;

void postprocess2d()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  int num_ps, i, j, k;
  double *pm, *pmstd;
  double temperature=1.0;
  double *lag, *Fcon_mean, *Fline2d_at_data_mean, *Fline2d_mean, *Trans2D_mean;
  void *posterior_sample, *post_model;
  double mean_lag, mean_lag_std, sum1, sum2;

  best_model_line2d = malloc(size_of_modeltype);
  best_model_std_line2d = malloc(size_of_modeltype);

// generate posterior sample
  temperature = 5.0;
  dnest_postprocess(temperature);

  if(thistask == roottask)
  {
    // initialize smooth work space
    smooth_init(n_vel_data);
    char fname[200];
    FILE *fp, *fcon, *fline, *ftran;

    // get number of lines in posterior sample file
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

    //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    //file for continuum reconstruction
    fcon = fopen("data/con_rec.txt", "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/con_rec.txt.\n");
      exit(0);
    }
    //file for line reconstruction
    fline = fopen("data/line2d_rec.txt", "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/line1d_rec.txt.\n");
      exit(0);
    }
    //file for transfer function
    ftran = fopen("data/tran2d_rec.txt", "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/tran_rec.txt.\n");
      exit(0);
    }

    // read number of lines
    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    lag = malloc(num_ps * sizeof(double));
    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);
    
    Fcon_mean = malloc(parset.n_con_recon*sizeof(double));
    Fline2d_at_data_mean = malloc(n_line_data * n_vel_data *sizeof(double));
    Fline2d_mean = malloc(parset.n_line_recon * parset.n_vel_recon *sizeof(double));
    Trans2D_mean = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));

    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;
    mean_lag = 0.0;
    for(j = 0; j<parset.n_con_recon; j++)
      Fcon_mean[j] = 0.0;
    for(j=0; j<n_line_data * n_vel_data; j++)
      Fline2d_at_data_mean[j] = 0.0;
    for(j=0; j<parset.n_line_recon * parset.n_vel_recon; j++)
      Fline2d_mean[j] = 0.0;
    for(j=0; j<parset.n_tau * parset.n_vel_recon; j++)
      Trans2D_mean[j] = 0.0;

    for(i=0; i<num_ps; i++)
    {
      // read lines
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      //store lines
      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);

      Fcon = Fcon_particles[which_particle_update];
    
      calculate_con_from_model(post_model + num_params_blr *sizeof(double));
      gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

      Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
      transfun_2d_cloud_direct(post_model, Vline_data, Trans2D_at_veldata, 
                                      n_vel_data, parset.flag_save_clouds);
      calculate_line2d_from_blrmodel(post_model, Tline_data, Vline_data, Trans2D_at_veldata, 
                                      Fline2d_at_data, n_line_data, n_vel_data);

      for(j=0; j<parset.n_con_recon; j++)
      {
        Fcon_mean[j] += Fcon[j];
      }
      for(j=0; j<n_line_data * n_vel_data; j++)
      {
        Fline2d_at_data_mean[j] += Fline2d_at_data[j];
      }
      
      sum1 = 0.0;
      sum2 = 0.0;
      for(j=0; j<parset.n_tau; j++)
      {
        for(k=0; k<n_vel_data; k++)
        {
          sum1 += Trans2D_at_veldata[j * n_vel_data + k] * TransTau[j];
          sum2 += Trans2D_at_veldata[j * n_vel_data + k];
        }
      }
      lag[i] = sum1/sum2;
      mean_lag += lag[i];

      if(gsl_rng_uniform(gsl_r) < 1.0)
      {
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%f %f\n", Tcon[j], Fcon[j]/con_scale);
        }
        fprintf(fcon, "\n");

        for(j=0; j<n_line_data; j++)
        {
          for(k=0; k<n_vel_data; k++)
          {
            fprintf(fline, "%f ", Fline2d_at_data[j * n_vel_data + k]/line_scale);
          }
          fprintf(fline, "\n");
        }
        fprintf(fline, "\n");

        for(j=0; j<n_line_data; j++)
        {
          for(k=0; k<n_vel_data; k++)
          {
            fprintf(ftran, "%f ", Trans2D_at_veldata[j * n_vel_data + k]);
          }
          fprintf(ftran, "\n");
        }
        fprintf(ftran, "\n");
      }
    }

    smooth_end();

    fclose(fp);
    fclose(fcon);
    fclose(fline);
    fclose(ftran);

    //
    smooth_init(parset.n_vel_recon);
    for(i=0; i<num_ps; i++)
    {
      memcpy(post_model, posterior_sample+i*size_of_modeltype, size_of_modeltype);

      Fcon = Fcon_particles[which_particle_update];
    
      calculate_con_from_model(post_model + num_params_blr *sizeof(double));
      gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

      transfun_2d_cloud_direct(post_model, TransV, Trans2D, parset.n_vel_recon, 0);
      calculate_line2d_from_blrmodel(post_model, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

      for(j=0; j<parset.n_line_recon * parset.n_vel_recon; j++)
      {
        Fline2d_mean[j] += Fline2d[j];
      }
      for(j=0; j<parset.n_tau * parset.n_vel_recon; j++)
      {
        Trans2D_mean[j] += Trans2D[j];
      }
    }
    smooth_end();

    for(j=0; j<parset.n_con_recon; j++)
      Fcon_mean[j] /= num_ps;
    for(j=0; j<n_line_data * n_vel_data; j++)
      Fline2d_at_data_mean[j] /= num_ps;
    for(j=0; j<parset.n_line_recon * parset.n_vel_recon; j++)
      Fline2d_mean[j] /= num_ps;
    for(j=0; j<parset.n_tau * parset.n_vel_recon; j++)
      Trans2D_mean[j] /= num_ps;
    
    sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tcon[i], Fcon_mean[i] / con_scale);
    }
    fclose(fp);

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
        fprintf(fp, "%f %f %f\n", Vline_data[j]*VelUnit, Tline_data[i],  Fline2d_at_data_mean[i*n_vel_data + j] / line_scale);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);

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
        fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, Tline[i],  Fline2d_mean[i*parset.n_vel_recon + j] / line_scale);
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
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<parset.n_vel_recon; j++)
      {
        fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, TransTau[i], Trans2D_mean[i*parset.n_vel_recon + j]);
      }
      fprintf(fp, "\n");
    }
    fclose(fp);

    mean_lag /= num_ps;
    mean_lag_std = 0.0;
    for(i=0; i<num_ps; i++)
    {
      mean_lag_std = (lag[i] - mean_lag) * (lag[i] - mean_lag);
    }
    if(num_ps > 1)
      mean_lag_std = sqrt(mean_lag_std/(num_ps -1.0));
    else
      mean_lag_std = 0.0;
    printf("Mean time lag: %f+-%f\n", mean_lag, mean_lag_std);

    pm = (double *)best_model_line2d;
    pmstd = (double *)best_model_std_line2d;
    for(j=0; j<num_params; j++)
    {
      pm[j] = pmstd[j] = 0.0;
    }
    for(i=0; i<num_ps; i++)
    {
      for(j =0; j<num_params; j++)
        pm[j] += *((double *)posterior_sample + i*num_params + j );
    }

    for(j=0; j<num_params; j++)
      pm[j] /= num_ps;

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
        pmstd[j] += pow( *((double *)posterior_sample + i*num_params + j ) - pm[j], 2.0 );
    }

    for(j=0; j<num_params; j++)
    {
      if(num_ps > 1)
        pmstd[j] = sqrt(pmstd[j]/(num_ps-1.0));
      else
        pmstd[j] = 0.0;
    }

    for(j = 0; j<num_params_blr + num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_line2d + j), *((double *)best_model_std_line2d+j) ); 
 
    free(lag);
    free(post_model);
    free(posterior_sample);
    free(Fcon_mean);
    free(Fline2d_at_data_mean);
    free(Fline2d_mean);
    free(Trans2D_mean);
  }
  return;
}

void reconstruct_line2d()
{
  char *argv[]={" "};

  reconstruct_line2d_init();
  
  smooth_init(n_vel_data);

  dnest_line2d(0, argv);

  smooth_end();

  postprocess2d();

  reconstruct_line2d_end();

  return;
}

void reconstruct_line2d_init()
{
  int i;
  double dT;
  /* set time array for continuum */
  Tcon_min = Tcon_data[0] - fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), parset.tau_max_set + (Tcon_data[0] - Tline_data[0]));
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), 10.0);
  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  TransTau = malloc(parset.n_tau * sizeof(double));
  //Trans2D_at_veldata = malloc(parset.n_tau * n_vel_data * sizeof(double));


  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans2D = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));
  
  Fline2d_at_data = malloc(n_line_data * n_vel_data * sizeof(double));

  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  
  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
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

  double vel_max_set = Vline_data[n_vel_data -1], vel_min_set = Vline_data[0];
  double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

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
  Trans2D_at_veldata_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans2D_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_data * sizeof(double));
    Trans2D_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_data * sizeof(double));
  }

  // only record gamma-distribution random number of clouds
  clouds_particles = malloc(parset.num_particles * sizeof(double *));
  clouds_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    clouds_particles[i] = malloc(parset.n_cloud_per_task * sizeof(double));
    clouds_particles_perturb[i] = malloc(parset.n_cloud_per_task * sizeof(double));
  }

  Fcon_particles = malloc(parset.num_particles * sizeof(double *));
  Fcon_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fcon_particles[i] = malloc(parset.n_con_recon * sizeof(double));
    Fcon_particles_perturb[i] = malloc(parset.n_con_recon * sizeof(double));
  }

  perturb_accept = malloc(parset.num_particles * sizeof(int));
  which_parameter_update_prev = malloc(parset.num_particles * sizeof(int));
  for(i=0; i<parset.num_particles; i++)
  {
    perturb_accept[i] = 0;
    which_parameter_update_prev[i] = -1;
  }

  prob_con_particles = malloc(parset.num_particles * sizeof(double));
  prob_con_particles_perturb = malloc(parset.num_particles * sizeof(double));
  return;
}

void reconstruct_line2d_end()
{
  free(Tline);
  free(Fline2d_at_data);
  free(Fline2d);

  free(TransTau);
  free(TransV);
  free(Trans2D);

  int i;
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans2D_at_veldata_particles[i]);
    free(Trans2D_at_veldata_particles_perturb[i]);
    free(clouds_particles[i]);
    free(clouds_particles_perturb[i]);
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(Trans2D_at_veldata_particles_perturb);
  free(clouds_particles);
  free(clouds_particles_perturb);
  free(Fcon_particles);
  free(Fcon_particles_perturb);
  
  free(perturb_accept);
  free(which_parameter_update_prev);
  free(prob_con_particles);
  free(prob_con_particles_perturb);

  free(par_fix);
  free(par_fix_val);

  free(best_model_line2d);
  free(best_model_std_line2d);
  return;
}

double prob_line2d(const void *model)
{
  double prob = 0.0, fcon, var2, dy;
  int i, param;
  double *pm = (double *)model;
  
  // if the previous perturb is accepted, store the previous Fcon at perturb stage;
  // otherwise, Fcon has no changes;
  // note that every time, only one parameter is updated, so that some 
  // perturb values keep unchanged. It is hard to track the updated parameters at last time
  if(perturb_accept[which_particle_update] == 1)
  {
    param = which_parameter_update_prev[which_particle_update];
    if(param >= num_params_blr)
    {
      memcpy(Fcon_particles[which_particle_update], Fcon_particles_perturb[which_particle_update], 
        parset.n_con_recon*sizeof(double));

      prob_con_particles[which_particle_update] = prob_con_particles_perturb[which_particle_update];
    }

    if(param < num_params_blr -1 )
    {
      memcpy(Trans2D_at_veldata_particles[which_particle_update], Trans2D_at_veldata_particles_perturb[which_particle_update], 
        parset.n_tau * n_vel_data * sizeof(double));
    }
    
  }

  if(which_parameter_update >= num_params_blr || which_parameter_update == -1)
  {

    Fcon = Fcon_particles_perturb[which_particle_update];
    calculate_con_from_model(model + num_params_blr*sizeof(double));
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

    for(i=0; i<n_con_data; i++)
    {
      fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
      var2 = Fcerrs_data[i] * Fcerrs_data[i] + exp(pm[num_params_blr]) * exp(num_params_blr);
      prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - 0.5* (log(2.0*PI) + log(var2) );
    }
    prob_con_particles_perturb[which_particle_update] = prob;
  }
  else
  {
    Fcon = Fcon_particles[which_particle_update];
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
    prob += prob_con_particles[which_particle_update];

    // update perturb values
    /*prob_con_particles_perturb[which_particle_update] = prob_con_particles[which_particle_update];
    memcpy(Fcon_particles_perturb[which_particle_update], Fcon_particles[which_particle_update], 
      parset.n_con_recon*sizeof(double));*/
  }
  
  // only update transfer function when BLR model is changed.
  // pm[11] only appears as errors
  if( (which_parameter_update < num_params_blr-1) || which_parameter_update == -1)
  {

    Trans2D_at_veldata = Trans2D_at_veldata_particles_perturb[which_particle_update];
    transfun_2d_cloud_direct(model, Vline_data, Trans2D_at_veldata, n_vel_data, 0);
    //memcpy(Trans2D_at_veldata_particles[which_particle_update], Trans2D_at_veldata, 
    //  n_vel_data*parset.n_tau * sizeof(double));
  }
  else
  {
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];

    //update perturb values
    /*memcpy(Trans2D_at_veldata_particles_perturb[which_particle_update], Trans2D_at_veldata_particles[which_particle_update], 
        parset.n_tau * n_vel_data * sizeof(double));*/
  }

  calculate_line2d_from_blrmodel(model, Tline_data, Vline_data, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data);

  for(i=0; i<n_line_data*n_vel_data; i++)
  {
    dy = Fline2d_data[i] - Fline2d_at_data[i] ;
    var2 = Flerrs2d_data[i]*Flerrs2d_data[i];
    var2 += exp(pm[num_params_blr-1])*exp(pm[num_params_blr-1]);
    prob += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
  }

  if(isnan(prob))
    prob = -DBL_MAX;

  if(which_parameter_update == -1)
  {
    memcpy(Fcon_particles[which_particle_update], Fcon_particles_perturb[which_particle_update], 
      parset.n_con_recon*sizeof(double));
    memcpy(Trans2D_at_veldata_particles[which_particle_update], Trans2D_at_veldata_particles_perturb[which_particle_update], 
        parset.n_tau * n_vel_data * sizeof(double));
    prob_con_particles[which_particle_update] = prob_con_particles_perturb[which_particle_update];
  }

  which_parameter_update_prev[which_particle_update] = which_parameter_update;
  return prob;
}

