/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
// header file for DNEST
#include "dnestvars.h"

#include "allvars.h"
#include "dnest_con.h"
#include "proto.h"

int num_params;

void *best_model_thismodel, *best_model_std_thismodel;

int dnest_con(int argc, char **argv)
{
  num_params = parset.n_con_recon + 3;
  size_of_modeltype = num_params * sizeof(double);
  best_model_thismodel = malloc(size_of_modeltype);
  best_model_std_thismodel = malloc(size_of_modeltype);

  /* setup functions used for dnest*/
  from_prior = from_prior_thismodel;
  log_likelihoods_cal = log_likelihoods_cal_thismodel;
  perturb = perturb_thismodel;
  print_particle = print_particle_thismodel;
  copy_model = copy_model_thismodel;
  create_model = create_model_thismodel;
  get_num_params = get_num_params_thismodel;
  copy_best_model = copy_best_model_thismodel;
  
  sprintf(options_file, "%s/%s", parset.file_dir, "src/OPTIONSCON");
  
  dnest(argc, argv);
  
  if(thistask == 0)
  {
    int j;
    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_thismodel + j), *((double *)best_model_std_thismodel+j)  ); 
  }
  
  return 0;
}

/*===========================================*/
// users responsible for following struct definitions

void from_prior_thismodel(const void *model)
{
  int i;
  double *pm = (double *)model;
  
  pm[0] = -3.0 + dnest_rand()*3.0;
  pm[1] = 2.0 + dnest_rand()*8.0;
  pm[2] = 0.0 + dnest_rand()*2.0;
  for(i=0; i<parset.n_con_recon; i++)
    pm[i+3] = dnest_randn();
}

double log_likelihoods_cal_thismodel(const void *model)
{
  double logL;
  logL = prob_con_variability(model);
  return logL;
}

double perturb_thismodel(const void *model)
{
  double *pm = (double *)model;
  double logH = 0.0;
  //printf("P1 %f %f\n", model->params[0], model->params[1]);
  int which = dnest_rand_int(num_params);
  if(which >= num_params || which < 0)
  {
    printf("# Error: Incorrect which.\n");
    exit(0);
  }

  switch(which)
  {
    case 0:
      pm[0] += 3.0*dnest_randh();
      wrap(&(pm[0]), -3.0, 0.0);
      break;
    
    case 1:
      pm[1] += 8.0*dnest_randh();
      wrap(&(pm[1]), 2.0, 10.0);
      break;

    case 2:
      pm[2] += 2.0*dnest_randh();
      wrap(&(pm[2]), 0.0, 2.0);
      break;

    default:
      logH -= (-0.5*pow(pm[which], 2.0) );
      pm[which] += 20.0*dnest_randh();
      wrap(&pm[which], -10.0, 10.0);
      logH += (-0.5*pow(pm[which], 2.0) );
      break;
  }
  //printf("P2 %f %f %d\n", model->params[0], model->params[1], which);
  return logH;
}

void print_particle_thismodel(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}

void copy_model_thismodel(const void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void *create_model_thismodel()
{
  return (void *)malloc( size_of_modeltype );
}

int get_num_params_thismodel()
{
  return num_params;
}

void copy_best_model_thismodel(const void *bm, const void *bm_std)
{
  memcpy(best_model_thismodel, bm, size_of_modeltype);
  memcpy(best_model_std_thismodel, bm_std, size_of_modeltype);
}