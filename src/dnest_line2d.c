/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_line2d.c
 *  \brief run dnest for 2d line analysis.
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
#include "dnest_line2d.h"
#include "proto.h"

/*!
 * this function does dnest samling.
 */
int dnest_line2d(int argc, char **argv)
{
  int i;

  num_params_blr = 12;
  num_params = parset.n_con_recon + num_params_blr + num_params_var;
  size_of_modeltype = num_params * sizeof(double);
  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  /* setup functions used for dnest*/
  from_prior = from_prior_line2d;
  log_likelihoods_cal = log_likelihoods_cal_line2d;
  log_likelihoods_cal_initial = log_likelihoods_cal_initial_line2d;
  perturb = perturb_line2d;
  print_particle = print_particle_line2d;
  copy_model = copy_model_line2d;
  create_model = create_model_line2d;
  get_num_params = get_num_params_line2d;
  
  set_par_range_model2d();
  /* setup the fixed parameters */
  set_par_fix(num_params_blr);
  for(i=num_params_blr; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  strcpy(options_file, dnest_options_file);

  if(parset.flag_postprc == 0)
  {
    force_update = 0;
    dnest(argc, argv);
  }
  
  return 0;
}

/*!
 * this function setups parameter ranges.
 */
void set_par_range_model2d()
{
  int i;

  /* setup parameter range, BLR parameters first */
  for(i=0; i<num_params_blr; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];
  }
  // variability parameters
  for(i=num_params_blr; i<num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[i-num_params_blr][0];
    par_range_model[i][1] = var_range_model[i-num_params_blr][1];
  }
  // continuum ligth curve values
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[num_params_var][0];
    par_range_model[i][1] = var_range_model[num_params_var][1];
  }

  return;
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_line2d(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params_blr + num_params_var; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }

  i=num_params_blr-1;
  pm[i] = par_range_model[i][0] - dnest_rand() * ( par_range_model[i][1] - par_range_model[0][0] )*0.01;
  i=num_params_blr;
  pm[i] = par_range_model[i][0] - dnest_rand() * ( par_range_model[i][1] - par_range_model[0][0] )*0.01;

  // cope with fixed parameters
  for(i=0; i<num_params_blr; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  
  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var+num_params_blr] = dnest_randn();
  
  which_parameter_update = -1;

  return;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_line2d(const void *model)
{
  double logL;
  logL = prob_line2d(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_initial_line2d(const void *model)
{
  double logL;
  logL = prob_initial_line2d(model);
  return logL;
}

/*!
 * this function perturbs parameters.
 */
double perturb_line2d(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which; 

  /* 
   * fixed parameters need not to update 
   * perturb important parameters more frequently
   */
  do
  {
    rnd = dnest_rand();
    if(rnd < 0.5)
      which = dnest_rand_int(num_params_blr + num_params_var);
    else
      which = dnest_rand_int(parset.n_con_recon) + num_params_blr + num_params_var;
  }while(par_fix[which] == 1);
 

  if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level_update = which_level_update > (size_levels-30)?(size_levels-30):which_level_update;
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

  if(which < num_params_blr)
  {
    // set an upper limit to move steps of systematic error parameters.
    if(which == num_params_blr-1 || which == num_params_blr )
       width = fmin(width, (par_range_model[which][1] - par_range_model[which][0])*0.01 );

    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr + num_params_var)
  {
    logH -= (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
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
 * this function prints out parameters.
 */
void print_particle_line2d(FILE *fp, const void *model)
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
 * this function copy model from src to dest.
 */
void copy_model_line2d(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
  return;
}

/*!
 * this function allocates a model memoery. 
 */
void *create_model_line2d()
{
  return (void *)malloc(size_of_modeltype);
}

/*!
 * this function returns number of parameters.
 */
int get_num_params_line2d()
{
  return num_params;
}
