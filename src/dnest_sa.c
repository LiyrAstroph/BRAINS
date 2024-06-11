/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*! \file dnest_sa.c
 *  \brief run dnest sampling for sa analysis.
 */

#ifdef SpecAstro

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
 #include <stddef.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#include "brains.h"

DNestFptrSet *fptrset_sa;

int dnest_sa(int argc, char **argv)
{
  int i;
  double logz_sa;
  char dnest_data_dir[BRAINS_MAX_STR_LENGTH];

  set_sa_blr_model();
  
  num_params_blr = 0;  /* RM parameter numbers */
  num_params_sa_blr = num_params_sa_blr_model + num_params_sa_extpar;
  num_params_blr_tot = num_params_sa_blr;
  num_params_sa = num_params_sa_blr;
  num_params = num_params_sa;

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

  fptrset_sa = dnest_malloc_fptrset();
  /* setup functions used for dnest*/
  fptrset_sa->from_prior = from_prior_sa;
  fptrset_sa->print_particle = print_particle_sa;
  fptrset_sa->restart_action = restart_action_sa;
  fptrset_sa->accept_action = accept_action_sa;
  fptrset_sa->kill_action = kill_action_sa;
  fptrset_sa->perturb = perturb_sa;
  fptrset_sa->read_particle = read_particle_sa;

  fptrset_sa->log_likelihoods_cal_initial = log_likelihoods_cal_initial_sa;
  fptrset_sa->log_likelihoods_cal_restart = log_likelihoods_cal_restart_sa;
  fptrset_sa->log_likelihoods_cal = log_likelihoods_cal_sa;

  set_par_range_sa();
  set_par_fix_sa_blrmodel();

  for(i=num_params_sa_blr_model; i<num_params; i++)
  {
    par_fix[i] = 0;
    par_fix_val[i] = -DBL_MAX;
  }

  /* fix DA */
  par_fix[num_params_sa_blr_model] = 1;
  par_fix_val[num_params_sa_blr_model] = log(550.0);

  /* fix FA */
  par_fix[num_params_sa_blr_model+2] = 1;
  par_fix_val[num_params_sa_blr_model+2] = log(1.0);
  
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
    
  print_par_names_sa();

  force_update = parset.flag_force_update;
  if(parset.flag_para_name != 1)
  {
    strcpy(dnest_data_dir,parset.file_dir);
    strcat(dnest_data_dir, "/data/");
    logz_sa = dnest(argc, argv, fptrset_sa, num_params, NULL, NULL, NULL, dnest_data_dir, dnest_options_file, NULL, NULL);
  }
  
  dnest_free_fptrset(fptrset_sa);

  return 0;
}

void set_par_range_sa()
{
  int i;

  /* setup parameter range, BLR parameters first */
  for(i=0; i<num_params_sa_blr_model; i++)
  {
    par_range_model[i][0] = sa_blr_range_model[i][0];
    par_range_model[i][1] = sa_blr_range_model[i][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }

  /* sa extra parameters */
  for(i=num_params_sa_blr_model; i<num_params_sa; i++)
  {
    par_range_model[i][0] = sa_extpar_range[i-num_params_sa_blr_model][0];
    par_range_model[i][1] = sa_extpar_range[i-num_params_sa_blr_model][1];

    par_prior_model[i] = UNIFORM;
    par_prior_gaussian[i][0] = 0.0;
    par_prior_gaussian[i][1] = 0.0;
  }
  return;
}

/*!
 *  print names and prior ranges for parameters 
 *
 */
void print_par_names_sa()
{
  if(thistask != roottask)
    return;

  int i, j;
  FILE *fp;
  char fname[BRAINS_MAX_STR_LENGTH], str_fmt[BRAINS_MAX_STR_LENGTH], str_name[BRAINS_MAX_STR_LENGTH];

  sprintf(fname, "%s/%s", parset.file_dir, "data/para_names_sa.txt");
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
  for(j=0; j<num_params_sa_blr_model; j++)
  {
    i++;
    strcpy(str_name, "\0");
    if(BLRmodel_sa_name!=NULL && BLRmodel_sa_name[i] != NULL)
    {
      strcpy(str_name, "SA_BLR_model_");
      strcat(str_name, BLRmodel_sa_name[i]);
    }
    else 
    {
      strcpy(str_name, "SA_BLR_model");
    }
    fprintf(fp, str_fmt, i, str_name, par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  
  for(j=0; j<num_params_sa_extpar; j++)
  {
    i++;
    strcpy(str_name, "SA_Ext_Par_");
    strcat(str_name, sa_extpar_name[j]);
    fprintf(fp, str_fmt, i, str_name, par_range_model[i][0], par_range_model[i][1], par_prior_model[i],
                            par_fix[i], par_fix_val[i], par_prior_gaussian[i][0], par_prior_gaussian[i][1]);
  }
  fclose(fp);
  return;
}

/*!
 * this function generates a sample from prior.
 */
void from_prior_sa(void *model)
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
double log_likelihoods_cal_initial_sa(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function calculates likelihood at initial step.
 */
double log_likelihoods_cal_restart_sa(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_sa(const void *model)
{
  double logL;
  logL = prob_sa(model);
  return logL;
}

/*!
 * this function prints out parameters.
 */
void print_particle_sa(FILE *fp, const void *model)
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
void read_particle_sa(FILE *fp, void *model)
{
  int j;
  double *psample = (double *)model;

  for(j=0; j < num_params; j++)
  {
    if(fscanf(fp, "%lf", psample+j) < 1)
    {
      printf("%f\n", *psample);
      fprintf(stderr, "#Error: Cannot read file sample_sa.txt.\n");
      exit(0);
    }
  }
  return;
}

/*!
 * this function perturbs parameters.
 */
double perturb_sa(void *model)
{
  double *pm = (double *)model;
  double logH = 0.0, limit1, limit2, width;
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
  which_level_update = dnest_get_which_level_update();
  size_levels = dnest_get_size_levels();
  which_level = which_level_update > (size_levels-10)?(size_levels-10):which_level_update;

  if( which_level > 0)
  {
    limit1 = limits[(which_level-1) * num_params *2 + which *2];
    limit2 = limits[(which_level-1) * num_params *2 + which *2 + 1];
    width = (limit2 - limit1);
  }
  else
  {
    limit1 = par_range_model[which][0];
    limit2 = par_range_model[which][1];
    width = (par_range_model[which][1] - par_range_model[which][0]);
  }
  width /= (2.35);

  if(par_prior_model[which] == GAUSSIAN)
  {
    logH -= (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
    pm[which] += dnest_randh() * width;
    dnest_wrap(&pm[which], par_range_model[which][0], par_range_model[which][1]);
    //dnest_wrap(&pm[which], limit1, limit2);
    logH += (-0.5*pow((pm[which] - par_prior_gaussian[which][0])/par_prior_gaussian[which][1], 2.0) );
  }
  else
  {
    pm[which] += dnest_randh() * width;
    dnest_wrap(&(pm[which]), par_range_model[which][0], par_range_model[which][1]);
    //dnest_wrap(&pm[which], limit1, limit2);
  }

  return logH;
}


/*!
 * this function calculates likelihood.
 */
double log_likelihoods_cal_sa_exam(const void *model)
{
  return 0.0;
}

/*
 * action when a particle's move is accepted.
 * usually store some values with no need to recompute.
 */
void accept_action_sa()
{  
  return;
}

/*
 * action when particle i is killed in cdnest sampling.
 * particle i_copy's properties is copyed to particle i. 
 */
void kill_action_sa(int i, int i_copy)
{
  return;
}

#endif