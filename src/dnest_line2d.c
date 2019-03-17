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


DNestFptrSet *fptrset_line2d;

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
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 2:
      num_params_blr_model = sizeof(BLRmodel2)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model2;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 3:
      num_params_blr_model = sizeof(BLRmodel3)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model3;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
    case 4:
      num_params_blr_model = sizeof(BLRmodel4)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model4;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      
      break;
    case 5:
      num_params_blr_model = sizeof(BLRmodel5)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model5;

      num_params_radial_samp=5;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 5;
      params_radial_samp[4] = 11; // black hole mass changes Schwarzschild radius

      break;
    
    case 6:
      num_params_blr_model = sizeof(BLRmodel6)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model6;

      num_params_radial_samp=4;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 10; // black hole mass changes Schwarzschild radius
      
      break;
    
    case 7:
      num_params_blr_model = sizeof(BLRmodel7)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model7;

      num_params_radial_samp=8;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;
      params_radial_samp[3] = 10;
      params_radial_samp[4] = 11;
      params_radial_samp[5] = 12;
      params_radial_samp[6] = 13;
      params_radial_samp[7] = 15; // black hole mass changes Schwarzschild radius
      
      break;

    case 8:
      num_params_blr_model = sizeof(BLRmodel8)/sizeof(double);
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model8;

      num_params_radial_samp=6;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 7;
      params_radial_samp[2] = 8;
      params_radial_samp[3] = 11;
      params_radial_samp[4] = 12;
      params_radial_samp[5] = 15;
      
      break;
    default:
      num_params_blr_model = 12;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model1;

      num_params_radial_samp=3;
      params_radial_samp = malloc(num_params_radial_samp * sizeof(int));
      params_radial_samp[0] = 2;
      params_radial_samp[1] = 3;
      params_radial_samp[2] = 4;

      break;
  }
  
  num_params_blr = num_params_blr_model + num_params_nlr_oiii + num_params_nlr + num_params_res + num_params_linecenter;
  num_params = parset.n_con_recon + num_params_blr + num_params_var;
  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  par_range_model = malloc( num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
    par_range_model[i] = malloc(2*sizeof(double));

  fptrset_line2d = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_line2d->from_prior = from_prior_line2d;
  fptrset_line2d->print_particle = print_particle_line2d;
  fptrset_line2d->restart_action = restart_action_2d;
  fptrset_line2d->accept_action = accept_action_2d;
  fptrset_line2d->kill_action = kill_action_2d;
  fptrset_line2d->perturb = perturb_line2d;
  fptrset_line2d->read_particle = read_particle_line2d;
  
  if(parset.flag_exam_prior != 1)
  {
    fptrset_line2d->log_likelihoods_cal_initial = log_likelihoods_cal_initial_line2d;
    fptrset_line2d->log_likelihoods_cal_restart = log_likelihoods_cal_restart_line2d;
    fptrset_line2d->log_likelihoods_cal = log_likelihoods_cal_line2d;
  }
  else
  {
    fptrset_line2d->log_likelihoods_cal_initial = log_likelihoods_cal_line2d_exam;
    fptrset_line2d->log_likelihoods_cal_restart = log_likelihoods_cal_line2d_exam;
    fptrset_line2d->log_likelihoods_cal = log_likelihoods_cal_line2d_exam;
  }
  
  set_par_range_model2d();

  /* setup the fixed parameters */
  set_par_fix(num_params_blr);

  /* setup the remaining paramters */
  for(i=num_params_blr; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }
  /* if flag_fixvar is true, fix continuum variation parameter sigma and tau */
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
  logz_line2d = dnest(argc, argv, fptrset_line2d, num_params, dnest_options_file);
  
  dnest_free_fptrset(fptrset_line2d);
  return 0;
}

/*!
 * this function setups parameter ranges.
 * 
 * The order of parameters is:
 *   I.   blr model.............(except systematic error)
 *   II.  narrow line...........(if flag_narrowline is true)
 *   III. spectral broadening...()
 *   IV.  line center...........(if flag_linecenter is true)
 *   V.   systematic error......()
 *   VI.  variability...........()
 *   VII. long-term trend.......()
 *   VIII.different trend.......()
 *   IX.  continuum light curve.()
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
  /* cope with narrow line of [OIII] doublet */
  for(i=num_params_blr_model-1; i<num_params_blr_model-1+num_params_nlr_oiii; i++)
  {
    par_range_model[i][0] = nlr_range_model[i - (num_params_blr_model-1)][0];
    par_range_model[i][1] = nlr_range_model[i - (num_params_blr_model-1)][1];
  }
  /* cope with narrow line */
  for(i=num_params_blr_model-1+num_params_nlr_oiii; i<num_params_blr-num_params_res-num_params_linecenter-1; i++)
  {
    par_range_model[i][0] = nlr_range_model[i - (num_params_blr_model+num_params_nlr_oiii-1)][0];
    par_range_model[i][1] = nlr_range_model[i - (num_params_blr_model+num_params_nlr_oiii-1)][1];
  }
  /* cope with spectral broadening */
  for(i=num_params_blr-num_params_res-num_params_linecenter-1; i<num_params_blr-num_params_linecenter-1; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;
  }
  /* cope with line center */
  for(i=num_params_blr-num_params_linecenter-1; i< num_params_blr-1; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;
  }
  /* the last is systematic error */
  i = num_params_blr-1;
  par_range_model[i][0] = blr_range_model[num_params_blr_model-1][0];
  par_range_model[i][1] = blr_range_model[num_params_blr_model-1][1];

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
  /* long-term trend of continuum */
  for(i=3 + num_params_blr; i< 4 + parset.flag_trend + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];
  }
  /* different trend in continuum and line */
  for(i=4 + parset.flag_trend + num_params_blr; i< num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[4 + i - (4 + parset.flag_trend + num_params_blr)][0];
    par_range_model[i][1] = var_range_model[4 + i - (4 + parset.flag_trend + num_params_blr)][1];
  }
  /* continuum ligth curve values */
  for(i=num_params_blr+num_params_var; i<num_params; i++)
  {
    par_range_model[i][0] = var_range_model[4+num_params_difftrend][0];
    par_range_model[i][1] = var_range_model[4+num_params_difftrend][1];
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
  /* cope with flux of narrow line and [OIII] doublet */
  for(i=num_params_blr_model-1; i<num_params_blr-num_params_res-num_params_linecenter-1-2; i++)
  {
    if(parset.flag_narrowline == 2) // Gaussian prior
    {
      pm[i] = dnest_randn();
      dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
    }
    else  // logrithmic prior
    {
      pm[i] = par_range_model[i][0] + dnest_rand() * ( par_range_model[i][1] - par_range_model[i][0]  );
    }
  }
  /* cope with width and shift of narrow line */
  for(i=num_params_blr_model-1+num_params_nlr_oiii+1; i<num_params_blr-num_params_res-num_params_linecenter-1; i++)
  {
    pm[i] = dnest_randn();
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  /* cope with spectral boradening */
  for(i=num_params_blr-num_params_res-num_params_linecenter-1; i<num_params_blr-num_params_linecenter-1; i++)
  {
    pm[i] = dnest_randn();
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  /* cope with line center */
  for(i=num_params_blr-num_params_linecenter-1; i< num_params_blr-1; i++)
  {
    pm[i] = dnest_randn();
    dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
  }
  /* systematic error */
  i=num_params_blr-1;
  pm[i] = par_range_model[i][1] - dnest_rand() * ( par_range_model[i][1] - par_range_model[0][0] );

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
  
  /* cope with fixed parameters */
  for(i=0; i<num_params_blr + num_params_var; i++)
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
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
  return;
}

/*!
 * This function read the particle from the file.
 */
void read_particle_line2d(FILE *fp, void *model)
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
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level_update = dnest_get_which_level_update();
  which_level = which_level_update > (size_levels-10)?(size_levels-10):which_level_update;

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
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr-num_params_res-num_params_linecenter-1-2) /* cope with flux of narrow line */
  {
    if(parset.flag_narrowline==2)  // Gaussian prior
    {
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += dnest_randh() * width;
      dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
      logH += (-0.5*pow(pm[which], 2.0) );
    }
    else  // logrithmic prior
    {
      pm[which] += dnest_randh() * width;
      dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    }
  }
  else if(which < num_params_blr-num_params_res-num_params_linecenter-1 )  /* cope with width and shift of narrow line */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_blr-1)  /* spectral broadening and line center */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_blr) /* systematic error */
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else if(which < num_params_blr + 3)  /* variability */
  {
    logH -= (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which]-var_param[which - num_params_blr])/var_param_std[which - num_params_blr], 2.0) );
  }
  else if(which < num_params_blr + 4 + parset.flag_trend)  /* long-term trend */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow(pm[which], 2.0) );
  }
  else if(which < num_params_blr + num_params_var)  /* different trend in continuum and line */
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
  }
  else  /* continuum light curve, Gaussian priors */
  {
    logH -= (-0.5*pow(pm[which], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
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

void accept_action_2d()
{
  int i, param, flag_cpy=0;
  double *ptemp;

  // the parameter previously updated
  param = which_parameter_update;
  // continuum parameter is updated
  if(param >= num_params_blr)
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
     * when this parameter is updated, Trans2D and Fline are unchanged.
     *      b) Fline is always changed, except for param = num_params_blr-1 or num_params_blr.
     */
    if(param < num_params_blr -1 )
    {
      ptemp = Trans2D_at_veldata_particles[which_particle_update];
      Trans2D_at_veldata_particles[which_particle_update] = Trans2D_at_veldata_particles_perturb[which_particle_update];
      Trans2D_at_veldata_particles_perturb[which_particle_update] = ptemp;

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

void kill_action_2d(int i, int i_copy)
{
  memcpy(Fcon_particles[i], Fcon_particles[i_copy], parset.n_con_recon * sizeof(double));
  memcpy(Fline_at_data_particles[i], Fline_at_data_particles[i_copy], n_line_data * sizeof(double));
  memcpy(con_q_particles[i], con_q_particles[i_copy], nq * sizeof(double));
  memcpy(Trans2D_at_veldata_particles[i], Trans2D_at_veldata_particles[i_copy], parset.n_tau * n_vel_data * sizeof(double));
  memcpy(clouds_particles[i], clouds_particles[i_copy], parset.n_cloud_per_task * sizeof(double)) ;
  return;
}
