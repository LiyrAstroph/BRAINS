/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*! \file dnest_sa1d.c
 *  \brief run dnest sampling for sa and 1d rm analysis.
 */

#ifdef SA

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#include "brains.h"

DNestFptrSet *fptrset_sa1d;

/* 
 * run denst for SA and 1D RM analysis.
 */
int dnest_sa1d(int argc, char **argv)
{
  int i;
  double logz_sa1d;
  
  if(parset.flag_sa_par_mutual == 0) /* SA and RM have the same BLR */
  {
    set_blr_model2d();
    num_params_sa_blr_model = 0;
  }
  else /* SA and RM have different BLRs but have the same mbh and inc. */
  {
    set_blr_model1d();
    set_sa_blr_model();
  }
  num_params_blr = num_params_blr_model + 2 + 1; /* include A, Ag, and line sys err */
  num_params_rm = parset.n_con_recon + num_params_var + num_params_blr;

  num_params_sa_blr = num_params_sa_blr_model + num_params_sa_extpar;
  num_params_sa = num_params_sa_blr;
  num_params = num_params_sa + num_params_rm;

  num_params_blr_tot = num_params_blr + num_params_sa_blr;

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  par_range_model = malloc( num_params * sizeof(double *));
  par_prior_gaussian = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
  {
    par_range_model[i] = malloc(2*sizeof(double));
    par_prior_gaussian[i] = malloc(2*sizeof(double));
  }
  par_prior_model = malloc( num_params * sizeof(int));

  fptrset_sa1d = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_sa1d->from_prior = from_prior_sa1d;
  fptrset_sa1d->print_particle = print_particle_sa1d;
  fptrset_sa1d->restart_action = restart_action_sa;
  fptrset_sa1d->accept_action = accept_action_sa1d;
  fptrset_sa1d->kill_action = kill_action_sa1d;
  fptrset_sa1d->perturb = perturb_sa1d;
  fptrset_sa1d->read_particle = read_particle_sa1d;

  fptrset_sa1d->log_likelihoods_cal_initial = log_likelihoods_cal_initial_sa1d;
  fptrset_sa1d->log_likelihoods_cal_restart = log_likelihoods_cal_restart_sa1d;
  fptrset_sa1d->log_likelihoods_cal = log_likelihoods_cal_sa1d;

  set_par_range_sa1d();
  print_par_names_sa1d();
  set_par_fix_blrmodel();

  exit(0);

  /* the rest parameters */
  for(i=num_params_blr_model; i<num_params; i++)
  {
    par_fix[i] = 0;
  }
  
  /* fix non-linear response */
  if(parset.flag_nonlinear !=1)
  {
    par_fix[num_params_blr-2] = 1;
    par_fix_val[num_params_blr-2] = 0.0;
  }

  /* fix systematic error of line */
  if(parset.flag_line_sys_err != 1)
  {
    par_fix[num_params_blr-1] = 1;
    par_fix_val[num_params_blr-1] = log(1.0);
  }

  /* fix mbh and inc of SA BLR if SA and RM have different BLR */
  if(parset.flag_sa_par_mutual != 0)
  {
    /* get the index for mbh and inclination parameters */
    set_idx_sa_par_mutual();
    /* mbh */
    par_fix[num_params_blr + idx_sa_par_mutual[0]] = 1;
    par_fix_val[num_params_blr + idx_sa_par_mutual[0]] = 0.0;
    
    /* inc */
    par_fix[num_params_blr + idx_sa_par_mutual[1]] = 1;
    par_fix_val[num_params_blr + idx_sa_par_mutual[1]] = 0.0;
  }

  /* fix FA */
  par_fix[num_params_blr_tot-1] = 1;
  par_fix_val[num_params_blr_tot-1] = log(1.0);

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err != 1)
  {
    par_fix[num_params_blr_tot] = 1;
    par_fix_val[num_params_blr_tot] = log(1.0);
  }

  /* fix continuum variation parameter sigma and tau if flag_fixvar is true */
  if(parset.flag_fixvar == 1)
  {
    par_fix[num_params_blr_tot + 1] = 1;
    par_fix_val[num_params_blr_tot + 1] = var_param[1];
    par_fix[num_params_blr_tot + 2] = 1;
    par_fix_val[num_params_blr_tot + 2] = var_param[2];
  }

  force_update = parset.flag_force_update;
  if(parset.flag_para_name != 1)
    logz_sa1d = dnest(argc, argv, fptrset_sa1d, num_params, dnest_options_file);
  
  dnest_free_fptrset(fptrset_sa1d);

  return 0;
}

void set_par_range_sa1d()
{
  int i;

  /* RM BLR parameters first */
  for(i=0; i<num_params_blr_model; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* RM response A and Ag, 2 parameters */ 
  for(i=num_params_blr_model; i<num_params_blr_model+2; i++)
  {
    par_range_model[i][0] = resp_range[i-num_params_blr_model][0];
    par_range_model[i][1] = resp_range[i-num_params_blr_model][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  /* systematic line error */
  i = num_params_blr -1;
  par_range_model[i][0] = sys_err_line_range[0];
  par_range_model[i][1] = sys_err_line_range[1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;


  /* SA BLR parameters */
  for(i=num_params_blr; i<num_params_blr + num_params_sa_blr_model; i++)
  {
    par_range_model[i][0] = sa_blr_range_model[i-num_params_blr][0];
    par_range_model[i][1] = sa_blr_range_model[i-num_params_blr][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  /* SA extra parameters */
  for(i=num_params_blr + num_params_sa_blr_model; i<num_params_blr_tot; i++)
  {
    par_range_model[i][0] = sa_extpar_range[i-(num_params_blr + num_params_sa_blr_model)][0];
    par_range_model[i][1] = sa_extpar_range[i-(num_params_blr + num_params_sa_blr_model)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* variability parameters */
  /* first systematic error */
  i = num_params_blr_tot;
  par_range_model[i][0] = var_range_model[0][0];
  par_range_model[i][1] = var_range_model[0][1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  /* drw parameters */
  for(i=num_params_blr_tot+1; i<num_params_drw + num_params_blr_tot; i++)
  {
    if(var_param_std[i-num_params_blr_tot] > 0.0)
    {
      par_range_model[i][0] = var_param[i-num_params_blr_tot] - 5.0 * var_param_std[i-num_params_blr_tot];
      par_range_model[i][1] = var_param[i-num_params_blr_tot] + 5.0 * var_param_std[i-num_params_blr_tot];

      /* make sure that the range lies within the initial range */
      par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[i-num_params_blr_tot][0]);
      par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[i-num_params_blr_tot][1]);

      par_prior_model[i] = GAUSSIAN;
      par_prior_gaussian[i][0] = var_param[i-num_params_blr_tot];
      par_prior_gaussian[i][1] = var_param_std[i-num_params_blr_tot];
    }
    else
    {
      par_range_model[i][0] = var_range_model[i-num_params_blr_tot][0];
      par_range_model[i][1] = var_range_model[i-num_params_blr_tot][1];

      par_prior_model[i] = UNIFORM;
      par_prior_gaussian[i][0] = 0.0;
      par_prior_gaussian[i][1] = 0.0;
    }
  }
  /* long-term trend */
  for(i=num_params_drw + num_params_blr_tot; i< num_params_drw + num_params_trend + num_params_blr_tot; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  /* different trend */
  for(i=num_params_drw + num_params_trend + num_params_blr_tot; i< num_params_var + num_params_blr_tot; i++)
  {
    par_range_model[i][0] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr_tot)][0];
    par_range_model[i][1] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr_tot)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* continuum light curve parameters */
  for(i=num_params_blr_tot+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[4+num_params_difftrend][0];
    par_range_model[i][1] = var_range_model[4+num_params_difftrend][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  return;
}

/*!
 *  print names and prior ranges for parameters 
 *
 */
void print_par_names_sa1d()
{
  if(thistask!= roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_sa1d.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  printf("# Print parameter name in %s\n", fname);

  i=-1;
  for(j=0; j<num_params_blr_model; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "BLR model", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "A", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "Ag", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sys_err_line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);

  for(j=0; j<num_params_sa_blr_model; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "SA BLR model", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<num_params_sa_extpar; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "SA Extra Par", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "sigmad", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  i++;
  fprintf(fp, "%4d %-15s %f %f %d\n", i, "taud", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  
  for(j=0; j<num_params_trend; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<num_params_difftrend; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "diff trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }

  for(j=0; j<parset.n_con_recon; j++)
  {
    i++;
    fprintf(fp, "%4d %-15s %f %f %d\n", i, "time series", par_range_model[i][0], par_range_model[i][1], par_prior_model[i]);
  }
  
  fclose(fp);
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_sa1d(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    if(par_prior_model[i] == GAUSSIAN)
    {
      pm[i] = dnest_randn()*par_prior_gaussian[i][1] + par_prior_gaussian[i][0];
      dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
    }
    else
    {
      pm[i] = par_range_model[i][0] + dnest_rand() * (par_range_model[i][1] - par_range_model[i][0]);
    }
  }

  /* cope with fixed parameters */
  for(i=0; i<num_params; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }

  which_parameter_update = -1;

  return;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_initial_sa1d(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_sa1d(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_sa1d(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function prints out parameters.
 */
void print_particle_sa1d(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

/*!
 * This function read the particle from the file.
 */
void read_particle_sa1d(FILE *fp, void *model)
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
 * this function perturbs parameters.
 */
double perturb_sa1d(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level, count_saves; 

  /* 
   * fixed parameters need not to update 
   * perturb important parameters more frequently
   */
  do
  {
    which = dnest_rand_int(num_params);
  }while(par_fix[which] == 1);
  
  which_parameter_update = which;
  
  /* level-dependent width */
  count_saves = dnest_get_count_saves();
  which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels-5)?(size_levels-5):which_level_update;

  if( which_level > 0 && count_saves > 1000)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = (limit2 - limit1);
  }
  else
  {
    limit1 = par_range_model[which][0];
    limit2 = par_range_model[which][1];
    width = (par_range_model[which][1] - par_range_model[which][0]);
  }
  width /= (5.0*2.35);

  if(par_prior_model[which] == GAUSSIAN)
  {
    logH -= (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
    pm[which] += dnest_randn() * width;
    //(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    dnest_wrap(&pm[which], limit1, limit2);
    logH += (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
  }
  else
  {
    pm[which] += dnest_randn() * width;
    //dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    dnest_wrap(&pm[which], limit1, limit2);
  }

  return logH;
}


/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_sa1d_exam(const void *model)
{
  return 0.0;
}

void accept_action_sa1d()
{
  /*int param;
  double *ptemp;

  // the parameter previously updated
  param = which_parameter_update;
  
  ptemp = phase_sa_particles[which_particle_update];
  phase_sa_particles[which_particle_update] = phase_sa_particles_perturb[which_particle_update];
  phase_sa_particles_perturb[which_particle_update] = ptemp;

  ptemp = Fline_sa_particles[which_particle_update];
  Fline_sa_particles[which_particle_update] = Fline_sa_particles_perturb[which_particle_update];
  Fline_sa_particles_perturb[which_particle_update] = ptemp;*/
  
  return;
}

void kill_action_sa1d(int i, int i_copy)
{
  memcpy(phase_sa_particles[i], phase_sa_particles[i_copy], n_vel_sa_data * n_base_sa_data * sizeof(double));
  memcpy(Fline_sa_particles[i], Fline_sa_particles[i_copy], n_vel_sa_data * n_epoch_sa_data * sizeof(double));
  return;
}

#endif
