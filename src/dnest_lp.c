/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_lp.c
 *  \brief run dnest for line profile fitting.
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


DNestFptrSet *fptrset_lp;

/*!
 * this function does dnest sampling.
 */
int dnest_lp(int argc, char **argv)
{
  int i, idx;
  char dnest_data_dir[BRAINS_MAX_STR_LENGTH];

  set_blr_model2d();

  num_params_blr = num_params_blr_model + num_params_nlr
                 + num_params_res + num_params_linecenter + 1; /* profile scale factor */
  num_params = num_params_blr;

  idx_linecenter = num_params_blr_model + num_params_nlr + num_params_res;

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

  fptrset_lp = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_lp->from_prior = from_prior_lp;
  fptrset_lp->print_particle = print_particle_lp;
  fptrset_lp->restart_action = restart_action_lp;
  fptrset_lp->accept_action = accept_action_lp;
  fptrset_lp->kill_action = kill_action_lp;
  fptrset_lp->perturb = perturb_lp;
  fptrset_lp->read_particle = read_particle_lp;

  if(parset.flag_exam_prior != 1)
  {
    fptrset_lp->log_likelihoods_cal_initial = log_likelihoods_cal_initial_lp;
    fptrset_lp->log_likelihoods_cal_restart = log_likelihoods_cal_restart_lp;
    fptrset_lp->log_likelihoods_cal = log_likelihoods_cal_lp;
  }
  else
  {
    fptrset_lp->log_likelihoods_cal_initial = log_likelihoods_cal_lp_exam;
    fptrset_lp->log_likelihoods_cal_restart = log_likelihoods_cal_lp_exam;
    fptrset_lp->log_likelihoods_cal = log_likelihoods_cal_lp_exam;
  }
  
  set_par_range_lp();
  /* setup the fixed parameters */
  set_par_fix_blrmodel();

  /* setup the remaining paramters */
  for(i=num_params_blr_model; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* fix mbh */
  idx = get_idx_mbh_from_blrmodel();
  par_fix[idx] = 1;
  par_fix_val[idx] =  log(10.0);
  
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

  /* cope with instrumental resolution */
  if(parset.flag_InstRes == 0)
  {
    par_fix[num_params_blr_model+num_params_nlr] = 1;
    par_fix_val[num_params_blr_model+num_params_nlr] = 0.0;
  }

  if(parset.flag_load_prior == 1)
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
  
  print_par_names_lp();

  force_update = parset.flag_force_update;
  if(parset.flag_para_name != 1)
  {
    strcpy(dnest_data_dir,parset.file_dir);
    strcat(dnest_data_dir, "/data/");
    logz_lp = dnest(argc, argv, fptrset_lp, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);
  }

  dnest_free_fptrset(fptrset_lp);
  return 0;
}

void set_par_range_lp()
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
  for(i=num_params_blr-num_params_linecenter-1; i< num_params_blr-1; i++)
  {
    par_range_model[i][0] = -10.0;
    par_range_model[i][1] =  10.0;

    par_prior_model[i] = GAUSSIAN;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 1.0;
  }

  /* the last is scale factor */
  i = num_params_blr-1;
  par_range_model[i][0] = -1.0;
  par_range_model[i][1] =  1.0;

  par_prior_model[i] = UNIFORM;
  par_prior_gaussian[i][0] = 0.0;
  par_prior_gaussian[i][1] = 0.0;

  return;
}

/*!
 *  print names and prior ranges for parameters 
 *
 */
void print_par_names_lp()
{
  if(thistask != roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH], str_name[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_lp.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  strcpy(str_fmt, "%4d %-28s %12.6f %12.6f %4d %4d %15.6e %15.6e %15.6e\n");

  printf("# Print parameter name in %s\n", fname);

  fprintf(fp, "#*************************************************\n");
  fprint_version(fp);
  fprintf(fp, "#*************************************************\n");

  fprintf(fp, "%4s %-28s %12s %12s %4s %4s %15s %15s %15s\n", "#", "Par", "Min", "Max", "Prior", "Fix", "Val", 
                                                              "Mean(Gau)", "Std(Gau)");

  i=-1;
  for(j=0; j<num_params_blr_model; j++)
  {
    i++;
    strcpy(str_name, "\0");
    if(BLRmodel_name!=NULL && BLRmodel_name[i] != NULL)
    {
      strcpy(str_name, "BLR_model_");
      strcat(str_name, BLRmodel_name[i]);
    }
    else 
    {
      strcpy(str_name, "BLR_model");
    }
    fprintf(fp, str_fmt, i, str_name, par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }

  for(j=0; j<num_params_nlr; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "narrow_line", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  for(j=0; j<num_params_res; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "line_broaden", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  for(j=0; j<num_params_linecenter; j++)
  {
    i++;
    fprintf(fp, str_fmt, i, "line_center", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }

  i++;
  fprintf(fp, str_fmt, i, "scale_factor", par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  
  fclose(fp);
  return;
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_lp(void *model)
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
  for(i=0; i<num_params; i++)
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
double log_likelihoods_cal_initial_lp(const void *model)
{
  double logL;
  logL = prob_lp(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_lp(const void *model)
{
  double logL;
  logL = prob_lp(model);
  return logL;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_lp(const void *model)
{
  double logL;
  logL = prob_lp(model);
  return logL;
}

/*!
 * this function prints out parameters.
 */
void print_particle_lp(FILE *fp, const void *model)
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
void read_particle_lp(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file samplelp.txt.\n");
      exit(0);
    }
  }
  return;
}

/*!
 * this function perturbs parameters.
 */
double perturb_lp(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width, rnd;
  int which, which_level, size_levels; 

  /* 
   * fixed parameters need not to update 
   * perturb important parameters more frequently
   */
  do
  {
    which = dnest_rand_int(num_params);
  }while(par_fix[which] == 1);
  
  which_parameter_update = which;
  
  /* level-dependent width */
 
  width = ( par_range_model[which][1] - par_range_model[which][0] );

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
double log_likelihoods_cal_lp_exam(const void *model)
{
  return 0.0;
}

void accept_action_lp()
{   
  return;
}

/*
 * action when particle i is killed in cdnest sampling.
 * particle i_copy's properties is copyed to particle i. 
 */
void kill_action_lp(int i, int i_copy)
{
  memcpy(Fline_at_data_particles[i], Fline_at_data_particles[i_copy], n_vel_data_ext  * sizeof(double));
  return;
}

void restart_action_lp(int iflag)
{
  return;
}