/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*! \file dnest_line1d.c
 *  \brief run dnest sampling for 1d line analysis.
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
 
// header file for DNEST
#include "dnestvars.h"

#include "allvars.h"
#include "dnest_line1d.h"
#include "proto.h"

DNestFptrSet *fptrset_line1d;

/*!
 * This function setup functions for dnest and run dnest.
 */
int dnest_line1d(int argc, char **argv)
{
  int i;
  
  switch(parset.flag_blrmodel)
  {
    case 1:
      num_params_blr_model = 9;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;


      break;
    case 2:
      num_params_blr_model = 9;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 3:
      num_params_blr_model = 9;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model3;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 4:
      num_params_blr_model = 9;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model3;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;

    case 5:
      num_params_blr_model = 12;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model5;

      num_params_radial_samp=4;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 5;
      break;

    case 6:
      num_params_blr_model = 11;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model6;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      break;

    case 7:
      num_params_blr_model = 16;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model7;

      num_params_radial_samp=7;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 10;
      params_radial_samp[4] = 11;
      params_radial_samp[5] = 12;
      params_radial_samp[6] = 13;
      break;

    default:
      num_params_blr_model = 9;
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      break;
  }
  
  num_params_blr = num_params_blr_model;
  num_params = parset.n_con_recon + num_params_var + num_params_blr;

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  fptrset_line1d = dnest_malloc_fptrset();

  /* setup functions used for dnest */
  fptrset_line1d->from_prior = from_prior_line1d;
  fptrset_line1d->print_particle = print_particle_line1d;
  fptrset_line1d->restart_action = restart_action_1d;
  fptrset_line1d->accept_action = accept_action_1d;
  fptrset_line1d->kill_action = kill_action_1d;
  fptrset_line1d->read_particle = read_particle_line1d;
  fptrset_line1d->perturb = perturb_line1d;

  if(parset.flag_exam_prior != 1)
  {
    fptrset_line1d->log_likelihoods_cal_initial = log_likelihoods_cal_initial_line1d;
    fptrset_line1d->log_likelihoods_cal_restart = log_likelihoods_cal_restart_line1d;
    fptrset_line1d->log_likelihoods_cal = log_likelihoods_cal_line1d;
  }
  else
  {
    fptrset_line1d->log_likelihoods_cal_initial = log_likelihoods_cal_line1d_exam;
    fptrset_line1d->log_likelihoods_cal_restart = log_likelihoods_cal_line1d_exam;
    fptrset_line1d->log_likelihoods_cal = log_likelihoods_cal_line1d_exam;
  }
  
  set_par_range_model1d();
  set_par_fix(num_params_blr);

  for(i=num_params_blr; i<num_params; i++)
    par_fix[i] = 0;

  /* fix continuum variation parameter sigma and tau if flag_fixvar is true */
  if(parset.flag_fixvar == 1)
  {
    par_fix[num_params_blr + 1] = 1;
    par_fix_val[num_params_blr + 1] = var_param[1];
    par_fix[num_params_blr + 2] = 1;
    par_fix_val[num_params_blr + 2] = var_param[2];
  }

  /* fix systematic error of line */
  if(parset.flag_line_sys_err != 1)
  {
    par_fix[num_params_blr-1] = 1;
    par_fix_val[num_params_blr-1] = log(1.0);
  }

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err != 1)
  {
    par_fix[num_params_blr] = 1;
    par_fix_val[num_params_blr] = log(1.0);
  }
  
  force_update = parset.flag_force_update;
  logz_line = dnest(argc, argv, fptrset_line1d, num_params, dnest_options_file);

  dnest_free_fptrset(fptrset_line1d);
  return 0;
}

/*!
 * this function set the parameter range.
 */
void set_par_range_model1d()
{
  int i;

  /* BLR parameters first */
  for(i=0; i<num_params_blr-1; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];
  }
  /* note that the last BLR parameters is the systematic error (1d) */
  i = num_params_blr -1;
  par_range_model[i][0] = blr_range_model[BLRmodel_size/sizeof(double)-1][0];
  par_range_model[i][1] = blr_range_model[BLRmodel_size/sizeof(double)-1][1];

  /* variability parameters */
  for(i=num_params_blr; i<3 + num_params_blr; i++)
  {
    if(var_param_std[i-num_params_blr] > 0.0)
    {
      par_range_model[i][0] = var_param[i-num_params_blr] - 5.0 * var_param_std[i-num_params_blr];
      par_range_model[i][1] = var_param[i-num_params_blr] + 5.0 * var_param_std[i-num_params_blr];

      /* make sure that the range lies within the initial range */
      par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[i-num_params_blr][0]);
      par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[i-num_params_blr][1]);
    }
    else
    {
      par_range_model[i][0] = var_range_model[i-num_params_blr][0];
      par_range_model[i][1] = var_range_model[i-num_params_blr][1];
    }
  }
  for(i=3 + num_params_blr; i< 4 + parset.flag_trend + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];
  }
  for(i=4 + parset.flag_trend + num_params_blr; i< num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[4 + i - (4 + parset.flag_trend + num_params_blr)][0];
    par_range_model[i][1] = var_range_model[4 + i - (4 + parset.flag_trend + num_params_blr)][1];
  }

  /* continuum light curve parameters */
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[5][0];
    par_range_model[i][1] = var_range_model[5][1];
  }
  return;
}

/*!
 * This function generates a sample from the prior.
 */
void from_prior_line1d(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_blr; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }

  /* variability parameters
   * use priors from continuum reconstruction.
   */
  for(i=num_params_blr; i<num_params_blr+3; i++)
  {
    pm[i] = dnest_randn()*var_param_std[i-num_params_blr] + var_param[i-num_params_blr];
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  /* long-term trend */
  for(i=num_params_blr+3; i<num_params_blr+ 4 + parset.flag_trend; i++)
  {
    pm[i] = dnest_randn();
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  /* different trend in continuum and line */
  for( i = num_params_blr+ 4 + parset.flag_trend; i< num_params_blr + num_params_var; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }
  
  /* continuum light curve */
  for(i=0; i<parset.n_con_recon; i++)
  {
    pm[i+num_params_var+num_params_blr] = dnest_randn();
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }

  /* cope with fixed parameters. */
  for(i=0; i<num_params_blr + num_params_var; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }

  /* all parameters need to update at the initial step */
  which_parameter_update = -1;
  return;
}

/*!
 * This function calculate log likelihood probability at the initial step.
 */
double log_likelihoods_cal_initial_line1d(const void *model)
{
  double logL;
  logL = prob_initial_line1d(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability at the initial step.
 */
double log_likelihoods_cal_restart_line1d(const void *model)
{
  double logL;
  logL = prob_restart_line1d(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability.
 */
double log_likelihoods_cal_line1d(const void *model)
{
  double logL;
  logL = prob_line1d(model);
  return logL;
}

/*!
 * This function print out the particle into the file.
 */
void print_particle_line1d(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

/*!
 * This function read the particle from the file.
 */
void read_particle_line1d(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < dnest_num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file %s.\n", options.sample_file);
      exit(0);
    }
  }
  return;
}

/*!
 * this function perturbs the parameters.
 */
double perturb_line1d(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level;

  /*
   * sample BLR and variability parameters more frequently; 
   * fixed parameter needs not to update.
   */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.5)
      which = dnest_rand_int(num_params_blr + num_params_var);
    else
      which = dnest_rand_int(parset.n_con_recon) + num_params_blr + num_params_var;

  }while(par_fix[which]==1);
  
  which_parameter_update = which;

  /* level-dependent width */
  which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels - 10)?(size_levels -10):which_level_update;

  if( which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }

  if(which < num_params_blr)  /* blr model */
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr + 3)  /* variability, Gaussion priors */
  {
    logH -= (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
  }
  else if(which < num_params_blr + 4 + parset.flag_trend) /* long-term trend */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_blr + num_params_var) /* different trend in continuum and line */
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else   /* continuum light curve, Gaussion prior */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  return logH;
}

/*!
 * This function calculate log likelihood probability.
 */
double log_likelihoods_cal_line1d_exam(const void *model)
{
  return 0.0;
}

void accept_action_1d()
{
  int param, i, flag_cpy=0;
  double *ptemp;

  // the parameter previously updated
  param = which_parameter_update;
  /* continuum parameter is updated */
  if( param >= num_params_blr )
  {
    /* 
     *note that (response) Fline is also changed as long as Fcon is changed.
     *num_params_blr-the parameter is the systematic error of continuum.
     *the change of this parameter also changes continuum reconstruction.
     */

    ptemp = Fcon_particles[which_particle_update];
    Fcon_particles[which_particle_update] = Fcon_particles_perturb[which_particle_update];
    Fcon_particles_perturb[which_particle_update] = ptemp;

    ptemp = Fline_at_data_particles[which_particle_update];
    Fline_at_data_particles[which_particle_update] = Fline_at_data_particles_perturb[which_particle_update];
    Fline_at_data_particles_perturb[which_particle_update] = ptemp;

    ptemp = con_q_particles[which_particle_update];
    con_q_particles[which_particle_update] = con_q_particles_perturb[which_particle_update];
    con_q_particles_perturb[which_particle_update] = ptemp;
  }
  else 
  {
    /* BLR parameter is updated 
     * Note a) that the (num_par_blr-1)-th parameter is systematic error of line.
     * when this parameter is updated, Trans1D and Fline are unchanged.
     *      b) Fline is always changed, except param = num_params_blr-1 or num_params_blr.
     */
    if( param < num_params_blr-1 )
    {

      ptemp = Trans1D_particles[which_particle_update];
      Trans1D_particles[which_particle_update] = Trans1D_particles_perturb[which_particle_update];
      Trans1D_particles_perturb[which_particle_update] = ptemp;

      ptemp = Fline_at_data_particles[which_particle_update];
      Fline_at_data_particles[which_particle_update] = Fline_at_data_particles_perturb[which_particle_update];
      Fline_at_data_particles_perturb[which_particle_update] = ptemp;

      /* when force_update is true, no need to store the perturbed value */
      if(force_update == 0)
      {
        for(i=0; i<num_params_radial_samp; i++)
        {
          if(param == params_radial_samp[i])
          {
            flag_cpy = 1;
            break;
          }
        }

        if(flag_cpy == 1)
        {
          ptemp = clouds_particles[which_particle_update];
          clouds_particles[which_particle_update] = clouds_particles_perturb[which_particle_update];
          clouds_particles_perturb[which_particle_update] = ptemp;
        }
      }
    }
  }  
  
  return;
}

void kill_action_1d(int i, int i_copy)
{
  memcpy(Fcon_particles[i], Fcon_particles[i_copy], parset.n_con_recon * sizeof(double));
  memcpy(Fline_at_data_particles[i], Fline_at_data_particles[i_copy], n_line_data * sizeof(double));
  memcpy(con_q_particles[i], con_q_particles[i_copy], nq * sizeof(double));
  memcpy(Trans1D_particles[i], Trans1D_particles[i_copy], parset.n_tau * sizeof(double));
  memcpy(clouds_particles[i], clouds_particles[i_copy], parset.n_cloud_per_task * sizeof(double)) ;
  return;
}
