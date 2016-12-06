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

/*!
 * This function setup functions for dnest and run dnest.
 */
int dnest_line1d(int argc, char **argv)
{
  int i;
  
  num_params_var = 4;
  num_params_blr = 9;
  num_params = parset.n_con_recon + num_params_var + num_params_blr;
  size_of_modeltype = num_params * sizeof(double);

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));

  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  /* setup functions used for dnest*/
  from_prior = from_prior_line1d;
  log_likelihoods_cal = log_likelihoods_cal_line1d;
  log_likelihoods_cal_initial = log_likelihoods_cal_initial_line1d;
  perturb = perturb_line1d;
  print_particle = print_particle_line1d;
  copy_model = copy_model_line1d;
  create_model = create_model_line1d;
  get_num_params = get_num_params_line1d;
  
  set_par_range_model1d();
  set_par_fix(num_params_blr);
  for(i=num_params_blr; i<num_params; i++)
    par_fix[i] = 0;

  strcpy(options_file, dnest_options_file);
  
  if(parset.flag_postprc == 0)
  {
    force_update = 0;
    dnest(argc, argv);
  }

  return 0;
}

void set_par_range_model1d()
{
  int i;

  for(i=0; i<num_params_blr-1; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];
  }
  i = num_params_blr -1;
  par_range_model[i][0] = blr_range_model[sizeof(BLRmodel)/sizeof(double)-1][0];
  par_range_model[i][1] = blr_range_model[sizeof(BLRmodel)/sizeof(double)-1][1];

  for(i=num_params_blr; i<num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[i-num_params_blr][0];
    par_range_model[i][1] = var_range_model[i-num_params_blr][1];
  }
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[num_params_var][0];
    par_range_model[i][1] = var_range_model[num_params_var][1];
  }
}

/*!
 * This function generate a sample from the prior.
 */
void from_prior_line1d(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_blr + num_params_var; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }

  i=num_params_blr-1;
  pm[i] = par_range_model[i][1] - dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0] )*0.01;
  i=num_params_blr;
  pm[i] = par_range_model[i][1] - dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0] )*0.01;

  for(i=0; i<num_params_blr; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  
  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var+num_params_blr] = dnest_randn();

  /* all parameters need to update at the initial step */
  which_parameter_update = -1;
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
 * This function calculate log likelihood probability at the initial step.
 */
double log_likelihoods_cal_initial_line1d(const void *model)
{
  double logL;
  logL = prob_initial_line1d(model);
  return logL;
}

double perturb_line1d(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which;

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

  if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }
  
  which_parameter_update = which;

  /* level-dependent width */
  which_level_update = which_level_update > (size_levels - 30)?(size_levels -30):which_level_update;
  which_level_update = which_level_update <0?0:which_level_update;

  if( which_level_update != 0)
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    width = ( par_range_model[which][1] - par_range_model[which][0] );
  }

  if(which < num_params_blr + num_params_var)
  {
    if(which == num_params_blr-1 || which == num_params_blr )
       width = fmin(width, (par_range_model[which][1] - par_range_model[which][0])*0.01 );

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
void print_particle_line1d(FILE *fp, const void *model)
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
void copy_model_line1d(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

/*!
 * This function create a model.
 */
void *create_model_line1d()
{
  return (void *)malloc( size_of_modeltype );
}

/*!
 * This function return the number of parameters.
 */
int get_num_params_line1d()
{
  return num_params;
}
