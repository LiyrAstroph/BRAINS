/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
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

#include "brains.h"


DNestFptrSet *fptrset_line2d;

/*!
 * this function does dnest sampling.
 */
int dnest_line2d(int argc, char **argv)
{
  int i;

  set_blr_model2d();
    
  num_params_blr = num_params_blr_model + num_params_nlr 
                 + num_params_res + num_params_linecenter + 2 + 1; /* include A, Ag, and line sys err */
  num_params_blr_tot = num_params_blr;
  num_params = parset.n_con_recon + num_params_blr + num_params_var;
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
  set_par_fix_blrmodel();

  /* setup the remaining paramters */
  for(i=num_params_blr_model; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* cope with narrow line */
  if(parset.flag_narrowline == 2)
  {
    /* flux */
    if(parset.flux_narrowline_err == 0.0)
    {
      par_fix[num_params_blr_model] = 1.0;
      par_fix_val[num_params_blr_model] = 0.0;
    }
  }
  if(parset.flag_narrowline >= 2)
  {
    /* width */
    if(parset.width_narrowline_err == 0.0)
    {
      par_fix[num_params_blr_model+1] = 1.0;
      par_fix_val[num_params_blr_model+1] = 0.0;
    } 
    /* shift */
    if(parset.shift_narrowline_err == 0.0)
    {
      par_fix[num_params_blr_model+2] = 1.0;
      par_fix_val[num_params_blr_model+2] = 0.0;
    }
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

  /* fix non-linear response */
  if(parset.flag_nonlinear !=1)
  {
    par_fix[num_params_blr-2] = 1;
    par_fix_val[num_params_blr-2] = 0.0;
  }
  
  print_par_names_model2d();

  force_update = parset.flag_force_update;
  if(parset.flag_para_name != 1)
    logz_line2d = dnest(argc, argv, fptrset_line2d, num_params, dnest_options_file);
  
  dnest_free_fptrset(fptrset_line2d);
  return 0;
}

/*!
 * this function setups parameter ranges.
 * 
 * The order of parameters is:                              \n
 *   I.   blr model.............()                          \n
 *   II.  narrow line...........(if flag_narrowline is true)\n
 *   III. spectral broadening...()                          \n
 *   IV.  line center...........(if flag_linecenter is true)\n
 *   V.   response A and Ag.....()                          \n
 *   VI.  systematic error......()                          \n
 *   VII. variability...........()                          \n
 *   VIII.long-term trend.......()                          \n
 *   IX.  different trend.......()                          \n
 *   X.   continuum light curve.()                          \n
 */
void set_par_range_model2d()
{
  int i;

  /* setup parameter range, BLR parameters first */
  for(i=0; i<num_params_blr_model; i++)
  {
    par_range_model[i][0] = blr_range_model[i][0];
    par_range_model[i][1] = blr_range_model[i][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* cope with narrow line */
  for(i=num_params_blr_model; i<num_params_blr_model + num_params_nlr; i++)
  {
    par_range_model[i][0] = nlr_range_model[i - num_params_blr_model][0];
    par_range_model[i][1] = nlr_range_model[i - num_params_blr_model][1];
    
    par_prior_model[i] = nlr_prior_model[i - num_params_blr_model];
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0; /* note that for logarithm prior of flux, this value is not used, so does not matter */
  }
  /* cope with spectral broadening */
  for(i=num_params_blr_model + num_params_nlr; i<num_params_blr_model + num_params_nlr + num_params_res; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  /* cope with line center */
  for(i=num_params_blr-num_params_linecenter-3; i< num_params_blr-3; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  
  for(i=num_params_blr-3; i<num_params_blr-1; i++)
  {
    par_range_model[i][0] = resp_range[i - (num_params_blr-3)][0];
    par_range_model[i][1] = resp_range[i - (num_params_blr-3)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* the last is systematic error */
  i = num_params_blr-1;
  par_range_model[i][0] = sys_err_line_range[0];
  par_range_model[i][1] = sys_err_line_range[1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  /* variability parameters */
   /* first systematic error */
  i = num_params_blr;
  par_range_model[i][0] = var_range_model[0][0];
  par_range_model[i][1] = var_range_model[0][1];

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  for(i=num_params_blr+1; i<num_params_drw + num_params_blr; i++)
  {
    if(var_param_std[i-num_params_blr] > 0.0)
    {
      par_range_model[i][0] = var_param[i-num_params_blr] - 5.0 * var_param_std[i-num_params_blr];
      par_range_model[i][1] = var_param[i-num_params_blr] + 5.0 * var_param_std[i-num_params_blr];

      /* make sure that the range lies within the initial range */
      par_range_model[i][0] = fmax(par_range_model[i][0], var_range_model[i-num_params_blr][0]);
      par_range_model[i][1] = fmin(par_range_model[i][1], var_range_model[i-num_params_blr][1]);

      par_prior_model[i] = GAUSSIAN;
      par_prior_gaussian[i][0] = var_param[i-num_params_blr];
      par_prior_gaussian[i][1] = var_param_std[i-num_params_blr];
    }
    else
    {
      par_range_model[i][0] = var_range_model[i-num_params_blr][0];
      par_range_model[i][1] = var_range_model[i-num_params_blr][1];

      par_prior_model[i] = UNIFORM;
      par_prior_gaussian[i][0] = 0.0;
      par_prior_gaussian[i][1] = 0.0;
    }
  }
  /* long-term trend of continuum */
  for(i=num_params_drw + num_params_blr; i< num_params_drw + num_params_trend + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }
  /* different trend in continuum and line */
  for(i=num_params_drw + num_params_trend + num_params_blr; i< num_params_var + num_params_blr; i++)
  {
    par_range_model[i][0] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr)][0];
    par_range_model[i][1] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_blr)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* continuum ligth curve values */
  for(i=num_params_blr+num_params_var; i<num_params; i++)
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
void print_par_names_model2d()
{
  if(thistask != roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_model2d.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  strcpy(str_fmt, "%4d %-15s %10.6f %10.6f %4d %4d %15.6e\n");

  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-15s %10s %10s %4s %4s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val");

  i=-1;
  for(j=0; j<num_params_blr_model; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "BLR model", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  for(j=0; j<num_params_nlr; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "narrow line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }
  for(j=0; j<num_params_res; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "line broaden", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }
  for(j=0; j<num_params_linecenter; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "line center", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }
  
  i++;
  fprintf(fp, str_fmt, i, "A", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);

  i++;
  fprintf(fp, str_fmt, i, "Ag", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);

  i++;
  fprintf(fp, str_fmt, i, "sys_err_line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);

  i++;
  fprintf(fp, str_fmt, i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  i++;
  fprintf(fp, str_fmt, i, "sigmad", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  i++;
  fprintf(fp, str_fmt, i, "taud", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  
  for(j=0; j<num_params_trend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  for(j=0; j<num_params_difftrend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "diff trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }

  for(j=0; j<parset.n_con_recon; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "time series", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i]);
  }
  
  fclose(fp);
}


/*!
 * this function generates a sample from prior.
 */
void from_prior_line2d(void *model)
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
    if(rnd < fmax(0.2, 1.0*(num_params_blr+num_params_var)/num_params))
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

  width /= 2.35;

  if(par_prior_model[which] == GAUSSIAN)
  {
    logH -= (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    logH += (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
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
  int param;
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

    ptemp = con_q_particles[which_particle_update];
    con_q_particles[which_particle_update] = con_q_particles_perturb[which_particle_update];
    con_q_particles_perturb[which_particle_update] = ptemp;

    if(force_update != 1)
    {
      ptemp = Fline_at_data_particles[which_particle_update];
      Fline_at_data_particles[which_particle_update] = Fline_at_data_particles_perturb[which_particle_update];
      Fline_at_data_particles_perturb[which_particle_update] = ptemp;
    }
  }
  else if( param < num_params_blr -1 && force_update != 1)
  {
    /* BLR parameter is updated 
     * Note a) that the (num_par_blr-1)-th parameter is systematic error of line.
     * when this parameter is updated, Trans2D and Fline are unchanged.
     *      b) Fline is always changed, except for param = num_params_blr-1.
     */
    
    {
      ptemp = TransTau_particles[which_particle_update];
      TransTau_particles[which_particle_update] = TransTau_particles_perturb[which_particle_update];
      TransTau_particles_perturb[which_particle_update] = ptemp;

      ptemp = Trans2D_at_veldata_particles[which_particle_update];
      Trans2D_at_veldata_particles[which_particle_update] = Trans2D_at_veldata_particles_perturb[which_particle_update];
      Trans2D_at_veldata_particles_perturb[which_particle_update] = ptemp;

      ptemp = Fline_at_data_particles[which_particle_update];
      Fline_at_data_particles[which_particle_update] = Fline_at_data_particles_perturb[which_particle_update];
      Fline_at_data_particles_perturb[which_particle_update] = ptemp;

    }
  } 
    
  return;
}

void kill_action_2d(int i, int i_copy)
{
  memcpy(Fcon_particles[i], Fcon_particles[i_copy], parset.n_con_recon * sizeof(double));
  memcpy(Fline_at_data_particles[i], Fline_at_data_particles[i_copy], n_line_data * sizeof(double));
  memcpy(con_q_particles[i], con_q_particles[i_copy], nq * sizeof(double));
  memcpy(TransTau_particles[i], TransTau_particles[i_copy], parset.n_tau*sizeof(double));
  memcpy(Trans2D_at_veldata_particles[i], Trans2D_at_veldata_particles[i_copy], parset.n_tau * n_vel_data * sizeof(double));
  return;
}
