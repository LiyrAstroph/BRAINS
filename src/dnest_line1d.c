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

  /* setup functions used for dnest*/
  from_prior = from_prior_line1d;
  log_likelihoods_cal = log_likelihoods_cal_line1d;
  log_likelihoods_cal_initial = log_likelihoods_cal_initial_line1d;
  perturb = perturb_line1d;
  print_particle = print_particle_line1d;
  copy_model = copy_model_line1d;
  create_model = create_model_line1d;
  get_num_params = get_num_params_line1d;
  
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


/*===========================================*/
// users responsible for following struct definitions
void from_prior_line1d(void *model)
{
  int i;
  double *pm = (double *)model;

  pm[0] = range_model[0].mu + dnest_rand()*( range_model[1].mu - range_model[0].mu );
  pm[1] = range_model[0].beta + dnest_rand()*( range_model[1].beta - range_model[0].beta );
  pm[2] = range_model[0].F + dnest_rand()*( range_model[1].F - range_model[0].F );
  pm[3] = range_model[0].inc + dnest_rand()*( range_model[1].inc - range_model[0].inc );
  pm[4] = range_model[0].opn + dnest_rand()*( range_model[1].opn - range_model[0].opn );
  pm[5] = range_model[0].A + dnest_rand()*( range_model[1].A - range_model[0].A );
  pm[6] = range_model[0].Ag + dnest_rand()*( range_model[1].Ag - range_model[0].Ag );
  pm[7] = range_model[0].k + dnest_rand()*( range_model[1].k - range_model[0].k );
  pm[8] = range_model[1].logse - dnest_rand()*( range_model[1].logse - range_model[0].logse )*0.01;

  for(i=0; i<num_params_blr; i++)
  {
    if(par_fix[i] == 1)
      pm[i] = par_fix_val[i];
  }
  
  pm[9] = var_range_model[0][1] - dnest_rand()*(var_range_model[0][1] - var_range_model[0][0]) * 0.01;
  pm[10] = var_range_model[1][0] + dnest_rand()* (var_range_model[1][1] - var_range_model[1][0]);
  pm[11] = var_range_model[2][0] + dnest_rand()* (var_range_model[2][1] - var_range_model[2][0]);
  pm[12] = var_range_model[3][0] + dnest_rand()* (var_range_model[3][1] - var_range_model[3][0]);
  for(i=0; i<parset.n_con_recon; i++)
    pm[i+num_params_var+num_params_blr] = dnest_randn();

  which_parameter_update = -1;

}

double log_likelihoods_cal_line1d(const void *model)
{
  double logL;
  logL = prob_line1d(model);
  return logL;
}

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

  which_level_update = which_level_update > 10?10:which_level_update;
  which_level_update = which_level_update <0?0:which_level_update;

  if( which_level_update != 0)
  {
    limit1 = limits[(which_level_update-1) * num_params *2 + which *2];
    limit2 = limits[(which_level_update-1) * num_params *2 + which *2 + 1];
    width = limit2 - limit1;
  }

  switch(which)
  {
  	case 0: // mu
      if(which_level_update == 0)
      {
        width = ( range_model[1].mu - range_model[0].mu );
      }

      pm[0] += dnest_randh() * width;
      wrap(&(pm[0]), range_model[0].mu, range_model[1].mu);
      break;
    
    case 1: // beta
      if(which_level_update == 0)
      {
        width =  ( range_model[1].beta - range_model[0].beta );
      }
      pm[1] += dnest_randh() * width;
      wrap(&(pm[1]), range_model[0].beta, range_model[1].beta);
      break;

    case 2: // F
      if(which_level_update == 0)
      {
        width = ( range_model[1].F - range_model[0].F );
      }
      pm[2] += dnest_randh() * width;
      wrap(&(pm[2]), range_model[0].F, range_model[1].F);
      break;

    case 3: // inclination
      if(which_level_update == 0)
      {
        width = ( range_model[1].inc - range_model[0].inc );
      }
      pm[3] += dnest_randh() * width;
      wrap(&(pm[3]), range_model[0].inc, range_model[1].inc);
      break;

    case 4: // openning angle
      if(which_level_update == 0)
      {
        width = ( range_model[1].opn - range_model[0].opn );
      }
      pm[4] += dnest_randh() * width;
      wrap(&(pm[4]), range_model[0].opn, range_model[1].opn);
      break;

    case 5: // A, response coefficient
      if(which_level_update == 0)
      {
        width = ( range_model[1].A - range_model[0].A );
      }
      pm[5] += dnest_randh() * width;
      wrap(&(pm[5]), range_model[0].A, range_model[1].A);
      break;

    case 6: // Ag, non-linearity
      if(which_level_update == 0)
      {
        width = ( range_model[1].Ag - range_model[0].Ag );
      }
      pm[6] += dnest_randh() * width;
      wrap(&(pm[6]), range_model[0].Ag, range_model[1].Ag);
      break;

    case 7: // k
      if(which_level_update == 0)
      {
        width = ( range_model[1].k - range_model[0].k );
      }
      pm[7] += dnest_randh() * width;
      wrap(&(pm[7]), range_model[0].k, range_model[1].k);
      break;

     case 8: // systematic error of line
      if(which_level_update == 0)
      {
        width =  ( range_model[1].logse - range_model[0].logse );
      }
      pm[8] += dnest_randh() * fmin(width, (range_model[1].logse - range_model[0].logse)*0.01 );
      wrap_limit(&(pm[which]), range_model[0].logse, range_model[1].logse);
      break;

    case 9: // systematic error of continuum
      if(which_level_update == 0)
      {
        width = var_range_model[0][1] - var_range_model[0][0];
      }
      pm[9] += dnest_randh() * fmin(width, (var_range_model[0][1] - var_range_model[0][0]) * 0.01);
      wrap(&(pm[9]), var_range_model[0][0], var_range_model[0][1]);
      break;
    
    case 10: // sigma
      if(which_level_update == 0)
      {
        width = var_range_model[1][1] - var_range_model[1][0];
      }
      pm[10] += dnest_randh() * width;
      wrap(&(pm[10]), var_range_model[1][0], var_range_model[1][1]);
      break;

    case 11: // tau
      if(which_level_update == 0)
      {
        width = var_range_model[2][1] - var_range_model[2][0];
      }
      pm[11] += dnest_randh() * width;
      wrap(&(pm[11]), var_range_model[2][0], var_range_model[2][1]);
      break;

    case 12: // mean value
      if(which_level_update == 0)
      {
        width = var_range_model[3][1] - var_range_model[3][0];
      }
      pm[12] += dnest_randh() * width;
      wrap(&(pm[12]), var_range_model[3][0], var_range_model[3][1]);
      break;

    default: // light curve points
      if(which_level_update == 0)
      {
        width = var_range_model[4][1] - var_range_model[4][0];;
      }
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += dnest_randh() * width;
      wrap(&pm[which], var_range_model[4][0], var_range_model[4][1]);
      logH += (-0.5*pow(pm[which], 2.0) );
      break;
  }
  return logH;
}

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

void copy_model_line1d(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void *create_model_line1d()
{
  return (void *)malloc( size_of_modeltype );
}

int get_num_params_line1d()
{
  return num_params;
}
