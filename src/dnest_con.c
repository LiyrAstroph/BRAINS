/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
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

#include "brains.h"

DNestFptrSet *fptrset_con;

/*!
 *  This function run dnest sampling for continuum.
 */
int dnest_con(int argc, char **argv)
{
  int i;
  char dnest_data_dir[BRAINS_MAX_STR_LENGTH];

  num_params = parset.n_con_recon + num_params_var;
  rnd_frac = fmax(0.2, 1.0*(num_params_var)/num_params);

  par_range_model = malloc( num_params * sizeof(double *));
  par_prior_gaussian = malloc(num_params * sizeof(double *));
  for(i=0; i<num_params; i++)
  {
    par_range_model[i] = malloc(2*sizeof(double));
    par_prior_gaussian[i] = malloc(2*sizeof(double));
  }
  par_prior_model = malloc( num_params * sizeof(int));
  

  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));
  
  fptrset_con = dnest_malloc_fptrset();

  /* setup functions used for dnest*/
  fptrset_con->from_prior = from_prior_con;
  fptrset_con->perturb = perturb_con;
  fptrset_con->print_particle = print_particle_con;
  fptrset_con->restart_action = restart_action_con;
  fptrset_con->accept_action = accept_action_con;
  fptrset_con->kill_action = kill_action_con;
  fptrset_con->read_particle = read_particle_con;
  
  if(parset.flag_exam_prior != 1)
  {
    fptrset_con->log_likelihoods_cal = log_likelihoods_cal_con;
    fptrset_con->log_likelihoods_cal_initial = log_likelihoods_cal_initial_con;
    fptrset_con->log_likelihoods_cal_restart = log_likelihoods_cal_restart_con;
  }
  else
  {
    fptrset_con->log_likelihoods_cal = log_likelihoods_cal_con_exam;
    fptrset_con->log_likelihoods_cal_initial = log_likelihoods_cal_con_exam;
    fptrset_con->log_likelihoods_cal_restart = log_likelihoods_cal_con_exam;
  }
  
  set_par_range_con();
  
  /* setup fixed parameters */
  for(i=0; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* fix systematic error of continuum */
  if(parset.flag_con_sys_err != 1)
  {
    par_fix[0] = 1;
    par_fix_val[0] = log(1.0);
  }
  
  /* if flag_dim != 0, no need to load priors */
  if(parset.flag_load_prior == 1 && parset.flag_dim == 0)
  {
    load_par_names(parset.prior_file);
    for(i=0; i<num_params; i++)
    {
      MPI_Bcast(par_range_model[i], 2, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
      MPI_Bcast(par_prior_gaussian[i], 2, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
    }
    MPI_Bcast(par_prior_model, num_params, MPI_INT, roottask, MPI_COMM_WORLD);
    MPI_Bcast(par_fix, num_params, MPI_INT, roottask, MPI_COMM_WORLD);
    MPI_Bcast(par_fix_val, num_params, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  }
  
  print_par_names_con();

  /* if not only print parameter name */
  if(parset.flag_para_name != 1)
  {
    strcpy(dnest_data_dir,parset.file_dir);
    strcat(dnest_data_dir, "/data/");
    logz_con = dnest(argc, argv, fptrset_con, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);
  }
  
  dnest_free_fptrset(fptrset_con);
  return 0;
}

/*!
 * this function set the parameter range.
 */
void set_par_range_con()
{
  int i;

  /* variability parameters */
  for(i=0; i<num_params_drw; i++)
  {
    par_range_model[i][0] = var_range_model[i][0];
    par_range_model[i][1] = var_range_model[i][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* parameters for long-term trend */
  for(i=num_params_drw; i<num_params_drw + num_params_trend; i++) 
  {
    par_range_model[i][0] = var_range_model[3][0];
    par_range_model[i][1] = var_range_model[3][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  
  /* paramter for response 
   * in case of 0D, num_params_resp = 0 
   * 1D and 2D call recontruct_con, in which resp parameters are needed */
  for(i=num_params_drw + num_params_trend; i<num_params_drw + num_params_trend + num_params_resp; i++)
  {
    par_range_model[i][0] = resp_range[i - (num_params_drw + num_params_trend)][0];
    par_range_model[i][1] = resp_range[i - (num_params_drw + num_params_trend)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* parameter for trend difference */
  for(i= num_params_drw + num_params_trend + num_params_resp; i<num_params_var; i++) 
  {
    par_range_model[i][0] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_resp)][0];
    par_range_model[i][1] = var_range_model[4 + i - (num_params_drw + num_params_trend + num_params_resp)][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* continuum light curve parameters */
  for(i=num_params_var; i<num_params; i++)
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
void print_par_names_con()
{
  if(thistask!= roottask)
    return;
  
  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_con.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  strcpy(str_fmt, "%4d %-28s %10.6f %10.6f %4d %4d %15.6e %15.6e %15.6e\n");

  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-28s %10s %10s %4s %4s %15s %15s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val", 
                                                              "Mean(Gau)", "Std(Gau)");
  i=0;
  fprintf(fp, str_fmt, i, "sys_err_con", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  i++;
  fprintf(fp, str_fmt, i, "sigmad", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  i++;
  fprintf(fp, str_fmt, i, "taud", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  
  for(j=0; j<num_params_trend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  
  /* in case of 0D, num_params_resp = 0 
   * 1D and 2D call recontruct_con, in which resp parameters are needed */
  for(j=0; j<num_params_resp; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "resp", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }

  for(j=0; j<num_params_difftrend; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "diff_trend", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }

  for(j=0; j<parset.n_con_recon; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "time_series", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  
  fclose(fp);
}

/*!
 * This function generate a sample from the prior.
 */
void from_prior_con(void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    if(par_prior_model[i] == GAUSSIAN )
    {
      pm[i] = dnest_randn() * par_prior_gaussian[i][1] + par_prior_gaussian[i][0];
      dnest_wrap(&pm[i], par_range_model[i][0], par_range_model[i][1]);
    }
    else
    {
      pm[i] = par_range_model[i][0] + dnest_rand()*(par_range_model[i][1] - par_range_model[i][0]);
    }
  }

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
  logL = prob_con_variability_semiseparable(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability at the initial step.
 */
double log_likelihoods_cal_initial_con(const void *model)
{
  double logL;
  logL = prob_con_variability_initial_semiseparable(model);
  return logL;
}

/*!
 * This function calculate log likelihood probability at the restart step.
 */
double log_likelihoods_cal_restart_con(const void *model)
{
  double logL;
  logL = prob_con_variability_initial_semiseparable(model);
  return logL;
}
/*!
 * This function generate a new move of parameters.
 */
double perturb_con(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level, size_levels;
  
  /* sample variability parameters more frequently */
  do
  {
    rnd = dnest_rand();
    if(rnd < rnd_frac)
      which = dnest_rand_int(num_params_var);
    else
      which = dnest_rand_int(parset.n_con_recon) + num_params_var;
  }while(par_fix[which] == 1);
  
  which_parameter_update = which;
  
  /* level-dependent width */
  which_level_update = dnest_get_which_level_update(); 
  size_levels = dnest_get_size_levels();
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
 * This function print the particle into the file.
 */
void print_particle_con(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%e ", pm[i] );
  }
  fprintf(fp, "\n");
}

/*!
 * This function read the particle from the file.
 */
void read_particle_con(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file sample.txt.\n");
      exit(0);
    }
  }
  return;
}

void restart_action_con(int iflag)
{

  return;
}

void accept_action_con()
{
  int param;

  param = which_parameter_update;
  /* only update prob when variability parameters are updated. */
  if(param < num_params_var)
  {
    prob_con_particles[which_particle_update] = prob_con_particles_perturb[which_particle_update];
  }
  return;
}

void kill_action_con(int i, int i_copy)
{
  prob_con_particles[i] = prob_con_particles[i_copy];
  return;
}


/*!
 * exam for prior
 */
double log_likelihoods_cal_con_exam(const void *model)
{
  return 0.0;
}