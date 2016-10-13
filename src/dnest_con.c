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


int dnest_con(int argc, char **argv)
{
  num_params_var = 4;
  num_params = parset.n_con_recon + num_params_var;
  size_of_modeltype = num_params * sizeof(double);
  
  /* setup functions used for dnest*/
  from_prior = from_prior_con;
  log_likelihoods_cal = log_likelihoods_cal_con;
  perturb = perturb_con;
  print_particle = print_particle_con;
  copy_model = copy_model_con;
  create_model = create_model_con;
  get_num_params = get_num_params_con;
  
  strcpy(options_file, dnest_options_file);
  
  dnest(argc, argv);
  
  return 0;
}

/*===========================================*/
// users responsible for following struct definitions

void from_prior_con(void *model)
{
  int i;
  double *pm = (double *)model;
  
  pm[0] = var_range_model[0][1] - dnest_rand()*(var_range_model[0][1] - var_range_model[0][0]) * 0.01;
  pm[1] = var_range_model[1][0] + dnest_rand()*(var_range_model[1][1] - var_range_model[1][0]);
  pm[2] = var_range_model[2][0] + dnest_rand()*(var_range_model[2][1] - var_range_model[2][0]);
  pm[3] = var_range_model[3][0] + dnest_rand()*(var_range_model[3][1] - var_range_model[3][0]);

  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var] = dnest_randn();
}

double log_likelihoods_cal_con(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}

double perturb_con(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width;
  
  int which = dnest_rand_int(num_params);
  if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }
  
  which_parameter_update = which;
  
  if(which_level_update !=0 || which_level_update < 10)
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }
  else
  {
    limit1 = limits[(10-1) * num_params *2 + which *2];
    limit2 = limits[(10-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }

  switch(which)
  {
    case 0:
      if(which_level_update == 0)
      {
        width = var_range_model[0][1] - var_range_model[0][0];
      }
      pm[0] += fmin(width, (var_range_model[0][1] - var_range_model[0][0]) * 0.01) * dnest_randh();
      wrap(&(pm[0]), var_range_model[0][0], var_range_model[0][1]);
      break;
    
    case 1:
      if(which_level_update == 0)
      {
        width = var_range_model[1][1] - var_range_model[1][0];
      }
      pm[1] += width*dnest_randh();
      wrap(&(pm[1]), var_range_model[1][0] , var_range_model[1][1] );
      break;

    case 2:
      if(which_level_update == 0)
      {
        width = var_range_model[2][1] - var_range_model[2][0];
      }
      pm[2] += width*dnest_randh();
      wrap(&(pm[2]), var_range_model[2][0], var_range_model[2][1]);
      break;

    case 3:
      if(which_level_update == 0)
      {
        width = var_range_model[3][1] - var_range_model[3][0];
      }
      pm[3] += width*dnest_randh();
      wrap(&(pm[3]), var_range_model[3][0], var_range_model[3][1]);
      break;

    default:
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

void copy_model_con(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void *create_model_con()
{
  return (void *)malloc( size_of_modeltype );
}

int get_num_params_con()
{
  return num_params;
}