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

  switch(parset.flag_blrmodel)
  {
    case 1:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 2:
      num_params_blr_model = sizeof(BLRmodel2)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model2;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 3:
      num_params_blr_model = sizeof(BLRmodel3)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model3;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 4:
      num_params_blr_model = sizeof(BLRmodel4)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model4;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      
      break;
    case 5:
      num_params_blr_model = sizeof(BLRmodel5)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model5;

      num_params_radial_samp=4;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 5;

      break;
    
    case 6:
      num_params_blr_model = sizeof(BLRmodel6)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model6;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      
      break;
    
    case 7:
      num_params_blr_model = sizeof(BLRmodel7)/sizeof(double);
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model7;

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
      num_params_blr_model = 12;
      transfun_2d_cloud_direct = transfun_2d_cloud_direct_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
  }
  
  num_params_blr = num_params_blr_model + num_params_nlr;
  num_params = parset.n_con_recon + num_params_blr + num_params_var;
  size_of_modeltype = num_params * sizeof(double);
  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  /* setup functions used for dnest*/
  from_prior = from_prior_line2d;
  print_particle = print_particle_line2d;
  get_num_params = get_num_params_line2d;
  restart_clouds = restart_clouds_2d;
  perturb = perturb_line2d;
  
  if(parset.flag_exam_prior != 1)
  {
    log_likelihoods_cal_initial = log_likelihoods_cal_initial_line2d;
    log_likelihoods_cal_restart = log_likelihoods_cal_restart_line2d;
    log_likelihoods_cal = log_likelihoods_cal_line2d;
  }
  else
  {
    log_likelihoods_cal_initial = log_likelihoods_cal_line2d_exam;
    log_likelihoods_cal_restart = log_likelihoods_cal_line2d_exam;
    log_likelihoods_cal = log_likelihoods_cal_line2d_exam;
  }
  
  set_par_range_model2d();
  /* setup the fixed parameters */
  set_par_fix(num_params_blr);
  for(i=num_params_blr; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }
  // fix continuum variation parameter sigma and tau
  if(parset.flag_fixvar == 1)
  {
    par_fix[num_params_blr + 1] = 1;
    par_fix_val[num_params_blr + 1] = var_param[1];
    par_fix[num_params_blr + 2] = 1;
    par_fix_val[num_params_blr + 2] = var_param[2];
  }

  strcpy(options_file, dnest_options_file);

  force_update = 0;
  dnest(argc, argv);
  
  return 0;
}

/*!
 * this function setups parameter ranges.
 */
void set_par_range_model2d()
{
  int i;

  /* setup parameter range, BLR parameters first */
  for(i=0; i<num_params_blr_model-1; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];
  }
  //cope with narrow line
  for(i=num_params_blr_model-1; i<num_params_blr-1; i++)
  {
    par_range_model[i][0] = nlr_range_model[i - (num_params_blr_model-1)][0];
    par_range_model[i][1] = nlr_range_model[i - (num_params_blr_model-1)][1];
  }
  // the last is systematic error
  i = num_params_blr-1;
  par_range_model[i][0] = blr_range_model[num_params_blr_model-1][0];
  par_range_model[i][1] = blr_range_model[num_params_blr_model-1][1];

  // variability parameters
  for(i=num_params_blr; i<3 + num_params_blr; i++)
  {
    //par_range_model[i][0] = var_range_model[i-num_params_blr][0];
    //par_range_model[i][1] = var_range_model[i-num_params_blr][1];
      par_range_model[i][0] = var_param[i-num_params_blr] - 5.0 * var_param_std[i-num_params_blr];
      par_range_model[i][1] = var_param[i-num_params_blr] + 5.0 * var_param_std[i-num_params_blr];
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

  // continuum ligth curve values
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[5][0];
    par_range_model[i][1] = var_range_model[5][1];
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

  for(i=0; i<num_params_blr_model -1 ; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }
  //cope with flux of narrow line
  for(i=num_params_blr_model-1; i<num_params_blr-1-2; i++)
  {
    if(parset.flag_narrowline == 2) // Gaussian prior
    {
      pm[i] = dnest_randn();
    }
    else  // logrithmic prior
    {
      pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
    }
  }
  // cope with width and shift of narrow line
  for(i=num_params_blr_model; i<num_params_blr-1; i++)
  {
    pm[i] = dnest_randn();
  }

  // systematic error
  i=num_params_blr-1;
  pm[i] = par_range_model[i][1] - dnest_rand() * ( par_range_model[i][1] - par_range_model[0][0] );

  // variability parameters
  // use priors from continuum reconstruction.
  for(i=num_params_blr; i<num_params_blr+3; i++)
  {
    pm[i] = dnest_randn()*var_param_std[i-num_params_blr] + var_param[i-num_params_blr];
    wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  for(i=num_params_blr+3; i<num_params_blr+ 4 + parset.flag_trend; i++)
  {
    pm[i] = dnest_randn();
  }
  for( i = num_params_blr+ 4 + parset.flag_trend; i< num_params_blr + num_params_var; i++)
  {
    pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
  }

  // cope with fixed parameters
  for(i=0; i<num_params_blr + num_params_var; i++)
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
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_initial_line2d(const void *model)
{
  double logL;
  logL = prob_initial_line2d(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_line2d(const void *model)
{
  double logL;
  logL = prob_restart_line2d(model);
  return logL;
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
 * this function returns number of parameters.
 */
int get_num_params_line2d()
{
  return num_params;
}

/*=======================================================================
 * model 1
 *=======================================================================
 */
/*!
 * this function perturbs parameters.
 */
double perturb_line2d(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level; 

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
 

  /*if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }*/
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level = which_level_update > (size_levels-50)?(size_levels-50):which_level_update;

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

  if(which < num_params_blr_model - 1)
  {
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr_model && parset.flag_narrowline > 1) // cope with flux of narrow line
  {
    if(parset.flag_narrowline==2)  // Gaussian prior
    {
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += dnest_randh() * width;
      wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
      logH += (-0.5*pow(pm[which], 2.0) );
    }
    else  // logrithmic prior
    {
      pm[which] += dnest_randh() * width;
      wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    }
  }
  else if(which < num_params_blr-1 && parset.flag_narrowline > 1)  // cope with width and shift of narrow line
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_blr) // systematic error
  {
    pm[which] += dnest_randh() * width;
    wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr + 4 + parset.flag_trend)
  {
    logH -= (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
    pm[which] += dnest_randh() * width;
    wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
  }
  else if(which < num_params_blr + num_params_var)
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
 * this function calculates likelihood.
 */
double log_likelihoods_cal_line2d_exam(const void *model)
{
  return 0.0;
}
