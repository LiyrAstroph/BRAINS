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
  num_params = parset.n_con_recon + num_params_var;
  size_of_modeltype = num_params * sizeof(double);
  
  /* setup functions used for dnest*/
  from_prior = from_prior_con;
  log_likelihoods_cal = log_likelihoods_cal_con;
  log_likelihoods_cal_initial = log_likelihoods_cal_initial_con;
  perturb = perturb_con;
  print_particle = print_particle_con;
  copy_model = copy_model_con;
  create_model = create_model_con;
  get_num_params = get_num_params_con;
  
  strcpy(options_file, dnest_options_file);
  
  if(parset.flag_postprc == 0)
  {
    dnest(argc, argv);
  }
  
  return 0;
}

/*!
 * This function generate a sample from the prior.
 */
void from_prior_con(void *model)
{
  int i;
  double *pm = (double *)model;
  
  pm[0] = var_range_model[0][1] - dnest_rand()*(var_range_model[0][1] - var_range_model[0][0]) * 0.01;
  pm[1] = var_range_model[1][0] + dnest_rand()*(var_range_model[1][1] - var_range_model[1][0]);
  pm[2] = var_range_model[2][0] + dnest_rand()*(var_range_model[2][1] - var_range_model[2][0]);
  pm[3] = dnest_randn(); //var_range_model[3][0] + dnest_rand()*(var_range_model[3][1] - var_range_model[3][0]);

  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var] = dnest_randn();

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
 * This function generate a new move of parameters.
 */
double perturb_con(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which;
  
  /* sample variability parameters more frequently */
  rnd = dnest_rand();
  if(rnd < 0.1)
    which = dnest_rand_int(num_params_var);
  else
    which = dnest_rand_int(parset.n_con_recon) + num_params_var;

  if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level_update = which_level_update > (size_levels - 30)?(size_levels-30):which_level_update;
  which_level_update = which_level_update <0?0:which_level_update;
  if( which_level_update != 0)
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }

  switch(which)
  {
    case 0: // systematic error of continuum; constrain the minimum step size 
      if(which_level_update == 0)
      {
        width = var_range_model[0][1] - var_range_model[0][0];
      }
      pm[which] += fmin(width, (var_range_model[0][1] - var_range_model[0][0]) * 0.1) * dnest_randh();
      wrap(&(pm[which]), var_range_model[0][0], var_range_model[0][1]);
      break;
    
    case 1: // sigma
      if(which_level_update == 0)
      {
        width = var_range_model[1][1] - var_range_model[1][0];
      }
      pm[which] += width*dnest_randh();
      wrap(&(pm[which]), var_range_model[1][0] , var_range_model[1][1] );
      break;

    case 2: // tau
      if(which_level_update == 0)
      {
        width = var_range_model[2][1] - var_range_model[2][0];
      }
      pm[which] += width*dnest_randh();
      wrap(&(pm[which]), var_range_model[2][0], var_range_model[2][1]);
      break;

    case 3: // mean value
      if(which_level_update == 0)
      {
        width = var_range_model[3][1] - var_range_model[3][0];
      }
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += width*dnest_randh();
      wrap(&(pm[which]), var_range_model[3][0], var_range_model[3][1]);
      logH += (-0.5*pow(pm[which], 2.0) );
      break;

    default: // light curve points
      if(which_level_update == 0)
      {
        width = var_range_model[4][1] - var_range_model[4][0];
      }
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += width*dnest_randh();
      wrap(&pm[which], var_range_model[4][0], var_range_model[4][1]);
      logH += (-0.5*pow(pm[which], 2.0) );
      break;
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
 * This function copy the model from src to dest.
 */
void copy_model_con(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

/*!
 * This function create a model.
 */
void *create_model_con()
{
  return (void *)malloc( size_of_modeltype );
}

/*!
 * This function return the number of parameters.
 */
int get_num_params_con()
{
  return num_params;
}
