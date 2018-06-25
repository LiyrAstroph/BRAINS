/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_con.c
 *  \brief run dnest for continuum analysis
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>
// header file for DNEST
#include "dnestvars.h"

#include "allvars.h"
#include "dnest_con.h"
#include "proto.h"

/*!
 *  This function run denst sampling for continuum.
 */
int dnest_con(int argc, char **argv)
{
  int i;

  num_params = parset.n_con_recon + num_params_var;
  size_of_modeltype = num_params * sizeof(double);

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  /* setup functions used for dnest*/
  from_prior = from_prior_con;
  perturb = perturb_con;
  print_particle = print_particle_con;
  get_num_params = get_num_params_con;
  restart_action = restart_action_con;
  
  if(parset.flag_exam_prior != 1)
  {
    log_likelihoods_cal = log_likelihoods_cal_con;
    log_likelihoods_cal_initial = log_likelihoods_cal_initial_con;
    log_likelihoods_cal_restart = log_likelihoods_cal_restart_con;
  }
  else
  {
    log_likelihoods_cal = log_likelihoods_cal_con_exam;
    log_likelihoods_cal_initial = log_likelihoods_cal_con_exam;
    log_likelihoods_cal_restart = log_likelihoods_cal_con_exam;
  }
  
  set_par_range_con();
  
  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
    par_fix[i] = 0;

  /* fix systematic error of continuum */
  par_fix[0] = 1;
  par_fix_val[0] = log(1.0);

  strcpy(options_file, dnest_options_file);
  dnest(argc, argv);
  
  return 0;
}

/*!
 * this function set the parameter range.
 */
void set_par_range_con()
{
  int i;

  /* variability parameters */
  for(i=0; i<3; i++)
  {
    par_range_model[i][0] = var_range_model[i][0];
    par_range_model[i][1] = var_range_model[i][1];
  }

  /* parameters for long-term trend */
  for(i=3; i<4 + parset.flag_trend; i++) 
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];
  }
  
  /* parameter for trend difference */
  for(i= 4 + parset.flag_trend; i<num_params_var; i++) 
  {
    par_range_model[i][0] = var_range_model[4 + i - (4 + parset.flag_trend)][0];
    par_range_model[i][1] = var_range_model[4 + i - (4 + parset.flag_trend)][1];
  }

  /* continuum light curve parameters */
  for(i=num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[5][0];
    par_range_model[i][1] = var_range_model[5][1];
  }
  return;
}

/*!
 * This function generate a sample from the prior.
 */
void from_prior_con(void *model)
{
  int i;
  double *pm = (double *)model;
  
  for(i=0; i<3; i++)
  {
    pm[i] = var_range_model[i][0] + dnest_rand()*(var_range_model[i][1] - var_range_model[i][0]);
  }
  for(i=3; i<4 + parset.flag_trend; i++)
  {
    pm[i] = dnest_randn();
  }
  for(i=4+parset.flag_trend; i<num_params_var; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
  }

  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var] = dnest_randn();

  for(i=0; i<num_params_var; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  /* all parameters need to update at the initial step */
  which_parameter_update = -1;
}

/*!
 * This function calculate log likelihood probability.
 */
double log_likelihoods_cal_con(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability at the initial step.
 */
double log_likelihoods_cal_initial_con(const void *model)
{
  double logL;
  logL = prob_con_variability_initial(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability at the restart step.
 */
double log_likelihoods_cal_restart_con(const void *model)
{
  double logL;
  logL = prob_con_variability_initial(model);
  return logL;
}
/*!
 * This function generate a new move of parameters.
 */
double perturb_con(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level;
  
  /* sample variability parameters more frequently */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.1)
      which = dnest_rand_int(num_params_var);
    else
      which = dnest_rand_int(parset.n_con_recon) + num_params_var;
  }while(par_fix[which] == 1);
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level = which_level_update > (size_levels - 10)?(size_levels-10):which_level_update;
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

  if(which < 3)
  {
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < 4 + parset.flag_trend)
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_var)
  {
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  return logH;
}

/*!
 * This function print the particle into the file.
 */
void print_particle_con(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}

/*!
 * This function return the number of parameters.
 */
int get_num_params_con()
{
  return num_params;
}

void restart_action_con(int iflag)
{
  return;
}


/*!
 * exam for prior
 */
double log_likelihoods_cal_con_exam(const void *model)
{
  return 0.0;
}