/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_line1d.c
 *  \brief reconstruct 1d line and recover BLR model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_interp.h>

#include "dnestvars.h"

#include "dnest_line1d.h"
#include "allvars.h"
#include "proto.h"


void *best_model_line1d, *best_model_std_line1d;

void postprocess1d()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  double temperature=1.0;
  double *pm, *pmstd;
  int num_ps, i, j;
  double *lag;
  void *posterior_sample, *post_model;
  double mean_lag, mean_lag_std, sum1, sum2;
  
  best_model_line1d = malloc(size_of_modeltype);
  best_model_std_line1d = malloc(size_of_modeltype);

// generate posterior sample 
  temperature = 1.0; 
  dnest_postprocess(temperature);

  if(thistask == roottask)
  {
    char fname[200];
    FILE *fp, *fcon, *fline, *ftran;
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
    sprintf(fname, "%s/%s", parset.file_dir, "data/line_rec.txt");
    fline = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran_rec.txt");
    ftran = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    lag = malloc(num_ps * sizeof(double));
    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);
    mean_lag = 0.0;
    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);
      
      Fcon = Fcon_particles[which_particle_update];
      calculate_con_from_model(post_model + num_params_blr *sizeof(double));
      gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

      Trans1D = Trans1D_particles[which_particle_update];
      transfun_1d_cloud_direct(post_model);
      calculate_line_from_blrmodel(post_model, Tline, Fline, parset.n_line_recon);

      sum1 = 0.0;
      sum2 = 0.0;
      for(j=0; j<parset.n_tau; j++)
      {
        sum1 += Trans1D[j] * TransTau[j];
        sum2 += Trans1D[j];
      }
      lag[i] = sum1/sum2;
      mean_lag += lag[i];

      if(i % (num_ps/10+1) == 0)
      {
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%f %f\n", Tcon[j], Fcon[j]/con_scale);
        }
        fprintf(fcon, "\n");

        for(j=0; j<parset.n_line_recon; j++)
        {
          fprintf(fline, "%f %f\n", Tline[j], Fline[j]/line_scale);
        }
        fprintf(fline, "\n");

        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran, "%f %f\n", TransTau[j], Trans1D[j]);
        }
        fprintf(ftran, "\n");
      }
    }
    fclose(fp);
    fclose(fcon);
    fclose(fline);
    fclose(ftran);

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

    pm = (double *)best_model_line1d;
    pmstd = (double *)best_model_std_line1d;
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

    double *pm = (double *)best_model_line1d;
    for(j = 0; j<num_params_blr+num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_line1d + j), 
                                             *((double *)best_model_std_line1d+j) );
    
    free(lag);
    free(post_model);
    free(posterior_sample);
  }
  return;
}

void reconstruct_line1d()
{
  char *argv[]={""};
  reconstruct_line1d_init();
// dnest run
  dnest_line1d(0, argv);
  postprocess1d();

  if(thistask == roottask)
  {
    FILE *fp;
    char fname[200];
    int i;
    double *pm = (double *)best_model_line1d;

    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;

    Fcon = Fcon_particles[which_particle_update];
    calculate_con_from_model(best_model_line1d + num_params_blr *sizeof(double));
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
    
    
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

    Trans1D = Trans1D_particles[which_particle_update];
    transfun_1d_cloud_direct(best_model_line1d);
    calculate_line_from_blrmodel(best_model_line1d, Tline, Fline, parset.n_line_recon);

    // output reconstructed line light curve
    sprintf(fname, "%s/%s", parset.file_dir, parset.pline_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_line_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tline[i], Fline[i] / line_scale);
    }
    fclose(fp);

    // output transfer function.
    sprintf(fname, "%s/%s", parset.file_dir, parset.tran_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_tau; i++)
    {
      fprintf(fp, "%f %f\n", TransTau[i], Trans1D[i]);
    }
    fclose(fp);
  }

  reconstruct_line1d_end();
}

void reconstruct_line1d_init()
{
  int i;
  double dT;

  /* set time array for continuum */
  Tcon_min = Tcon_data[0] - fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), parset.tau_max_set + (Tcon_data[0] - Tline_data[0]));
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), 10.0);
  Tcon_max = fmax(Tcon_max, Tline_data[n_line_data -1]);  /* The time span shoud cover that of the emission line data */

  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  TransTau = malloc(parset.n_tau * sizeof(double));
  //Trans1D = malloc(parset.n_tau * sizeof(double));

  //Fline_at_data = malloc(n_line_data * sizeof(double));

  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Flerrs = malloc(parset.n_line_recon * sizeof(double));

  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  Tline_min = fmax(Tline_min, Tcon_min + parset.tau_max_set);
  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  Tline_max = fmin(Tline_max, Tcon_max);  /* The time span should be within that of the continuum */
  
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

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS1D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  Trans1D_particles = malloc(parset.num_particles * sizeof(double *));
  Trans1D_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans1D_particles[i] = malloc(parset.n_tau * sizeof(double));
    Trans1D_particles_perturb[i] = malloc(parset.n_tau * sizeof(double));
  }

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_line_data * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_line_data * sizeof(double));
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
  prob_line_particles = malloc(parset.num_particles * sizeof(double));
  prob_line_particles_perturb = malloc(parset.num_particles * sizeof(double));
}

void reconstruct_line1d_end()
{
  free(TransTau);
  //free(Fline_at_data);

  free(Tline);
  free(Fline);
  free(Flerrs);

  int i;
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans1D_particles[i]);
    free(Trans1D_particles_perturb[i]);
    free(clouds_particles[i]);
    free(clouds_particles_perturb[i]);
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
  }
  free(Trans1D_particles);
  free(Trans1D_particles_perturb);
  free(clouds_particles);
  free(clouds_particles_perturb);
  free(Fcon_particles);
  free(Fcon_particles_perturb);

  free(perturb_accept);
  free(which_parameter_update_prev);
  free(prob_con_particles);
  free(prob_con_particles_perturb);

  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);
  free(prob_line_particles);
  free(prob_line_particles_perturb);

  free(par_fix);
  free(par_fix_val);
  free(best_model_line1d);
  free(best_model_std_line1d);

  if(thistask == roottask)
  {
    printf("Ends reconstruct_line1d.\n");
  }
}


double prob_line1d(const void *model)
{
  double prob = 0.0, prob_line=0.0, fcon, var2, dy;
  int i, param;
  double *pm = (double *)model;
  
  // if the previous perturb is accepted, store the previous perturb values, otherwise, no changes;
  if(perturb_accept[which_particle_update] == 1)
  { 
    param = which_parameter_update_prev[which_particle_update];

    /* continuum parameter is updated */
    if( param >= num_params_blr )
    {
      /* probability is always changed. */
      prob_con_particles[which_particle_update] = prob_con_particles_perturb[which_particle_update];

      /* the num_params_blr-th parameter is systematic error of continuum, which 
       * only appear at the stage of calculating likelihood probability.
       * when this paramete is updated, Fcon is unchanged.  
       */
      if( param > num_params_blr )
        memcpy(Fcon_particles[which_particle_update], Fcon_particles_perturb[which_particle_update], 
          parset.n_con_recon*sizeof(double));
    }
    else 
    {
      /* BLR parameter is udated 
     * Note that the (num_par_blr-1)-th parameter is systematic error of line.
     * when this parameter is updated, Trans1D and Fline are unchanged.
     */
      prob_line_particles[which_particle_update] = prob_line_particles_perturb[which_particle_update];
      if( param < num_params_blr-1 )
      {
        memcpy(Trans1D_particles[which_particle_update], Trans1D_particles_perturb[which_particle_update], 
          parset.n_tau * sizeof(double));

        memcpy(Fline_at_data_particles[which_particle_update], Fline_at_data_particles_perturb[which_particle_update], 
          n_line_data * sizeof(double));
      }
    }  
  }

  /* only update continuum reconstruction when the corresponding parameters are updated
   * or force to update (which_parameter_update = -1) 
   */
  if((which_parameter_update >= num_params_blr ) || which_parameter_update == -1)
  {
    /* the num_params_blr-th parameter is systematic error of continuum, which 
     * only appears at the stage of calculating likelihood probability.
     * when this paramete is updated, no need to re-calculate the contionuum.  
     */
    if( which_parameter_update > num_params_blr || which_parameter_update == -1 )
    {
      Fcon = Fcon_particles_perturb[which_particle_update];
      calculate_con_from_model(model + num_params_blr*sizeof(double));
    }
    else  /* if systematic error of continuum is updated */
    {
      Fcon = Fcon_particles[which_particle_update];
    }
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

    for(i=0; i<n_con_data; i++)
    {
      fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
      var2 = Fcerrs_data[i] * Fcerrs_data[i] + exp(pm[num_params_blr]) * exp(num_params_blr);
      prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - 0.5* (log(2.0*PI) + log(var2) );
    }

    prob_con_particles_perturb[which_particle_update] = prob;

  }
  else /* continuum has no change, use the previous values */
  {
    Fcon = Fcon_particles[which_particle_update];
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
    prob += prob_con_particles[which_particle_update];
  }
  
  /* only update transfer function when BLR model is changed
   * or forced to update (which_parameter_update = -1)
   * Trans1D is a pointer to the transfer function
   */
  if( (which_parameter_update < num_params_blr-1) || which_parameter_update == -1)
  {
    /* re-point */
    Trans1D = Trans1D_particles_perturb[which_particle_update]; 
    transfun_1d_cloud_direct(model);
  }
  else
  {
    /* re-point */
    Trans1D = Trans1D_particles[which_particle_update];
  }

  /* no need to calculate line when only systematic error parameters are updated.
   * otherwise, always need to calculate line.
   */
  if( (which_parameter_update!= num_params_blr-1 && which_parameter_update!= num_params_blr) || which_parameter_update == -1 )
  {
    /* re-point */
    Fline_at_data = Fline_at_data_particles_perturb[which_particle_update];
    calculate_line_from_blrmodel(model, Tline_data, Fline_at_data, n_line_data);

    for(i=0; i<n_line_data; i++)
    {
      dy = Fline_data[i] - Fline_at_data[i] ;
      var2 = Flerrs_data[i]*Flerrs_data[i];
      var2 += exp(pm[num_params_blr-1]) * exp(pm[num_params_blr-1]);
      prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
    }
    /* backup prob_line */
    prob_line_particles_perturb[which_particle_update] = prob_line;
  }
  else
  {
    /* continuum systematic error is updated */
    if( which_parameter_update == num_params_blr )
    {
      prob_line = prob_line_particles[which_particle_update];
    }
    else  /* line systematic error is updated */
    {
      /* re-point */
      Fline_at_data = Fline_at_data_particles[which_particle_update];
      for(i=0; i<n_line_data; i++)
      {
        dy = Fline_data[i] - Fline_at_data[i] ;
        var2 = Flerrs_data[i]*Flerrs_data[i];
        var2 += exp(pm[num_params_blr-1]) * exp(pm[num_params_blr-1]);
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
      /* backup prob_line */
      prob_line_particles_perturb[which_particle_update] = prob_line;
    }
  }

  // force to update
  if(which_parameter_update == -1)
  {
    memcpy(Fcon_particles[which_particle_update], Fcon_particles_perturb[which_particle_update], 
      parset.n_con_recon*sizeof(double));

    memcpy(Trans1D_particles[which_particle_update], Trans1D_particles_perturb[which_particle_update], 
      parset.n_tau * sizeof(double));

    memcpy(Fline_at_data_particles[which_particle_update], Fline_at_data_particles_perturb[which_particle_update],
      n_line_data * sizeof(double));

    prob_con_particles[which_particle_update] = prob_con_particles_perturb[which_particle_update];
    prob_line_particles[which_particle_update] = prob_line_particles_perturb[which_particle_update];
  }

  which_parameter_update_prev[which_particle_update] = which_parameter_update;

  /* add up line probability */
  prob += prob_line;
  return prob;
}