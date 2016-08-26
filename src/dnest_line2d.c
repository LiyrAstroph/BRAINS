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
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

// header file for DNEST
#include "dnestvars.h"

#include "allvars.h"
#include "dnest_line2d.h"
#include "proto.h"

int num_params;

void *best_model_line2d, *best_model_std_line2d;

int dnest_line2d(int argc, char **argv)
{
  num_params = parset.n_con_recon + 3 + 11;
  size_of_modeltype = num_params * sizeof(double);
  best_model_line2d = malloc(size_of_modeltype);
  best_model_std_line2d = malloc(size_of_modeltype);

  /* setup functions used for dnest*/
  from_prior = from_prior_line2d;
  log_likelihoods_cal = log_likelihoods_cal_line2d;
  perturb = perturb_line2d;
  print_particle = print_particle_line2d;
  copy_model = copy_model_line2d;
  create_model = create_model_line2d;
  get_num_params = get_num_params_line2d;
  copy_best_model = copy_best_model_line2d;
  
  strcpy(options_file, dnest_options_file);
  dnest(argc, argv);
  
  if(thistask == 0)
  {
    int j;
    for(j = 0; j<num_params; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_line2d + j), *((double *)best_model_std_line2d+j) ); 
  }

  return 0;

}

/*===========================================*/
// users responsible for following struct definitions

void from_prior_line2d(void *model)
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
  pm[8] = range_model[0].mbh + dnest_rand()*( range_model[1].mbh - range_model[0].mbh );
  pm[9] = range_model[0].lambda + dnest_rand()*( range_model[1].lambda - range_model[0].lambda );
  pm[10] = range_model[0].q + dnest_rand()*( range_model[1].q - range_model[0].q );
  
  pm[11] = -3.0 + dnest_rand()*3.0;
  pm[12] =  2.0 + dnest_rand()*8.0;
  pm[13] = 0.0 + dnest_rand()*2.0;
  for(i=0; i<parset.n_con_recon; i++)
    pm[i+3+11] = dnest_randn();
  
  which_parameter_update = -1; // force to update the clouds's ridal distribution
}

double log_likelihoods_cal_line2d(const void *model)
{
  double logL;
  logL = prob_line2d(model);
  return logL;
}

double perturb_line2d(void *model)
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
  
  which_parameter_update = which;
  
  switch(which)
  {
  	case 0:
      pm[0] += dnest_randh() * ( range_model[1].mu - range_model[0].mu );
      wrap(&(pm[0]), range_model[0].mu, range_model[1].mu);
      break;
    
    case 1:
      pm[1] += dnest_randh() * ( range_model[1].beta - range_model[0].beta );
      wrap(&(pm[1]), range_model[0].beta, range_model[1].beta);
      break;

    case 2:
      pm[2] += dnest_randh() * ( range_model[1].F - range_model[0].F );
      wrap(&(pm[2]), range_model[0].F, range_model[1].F);
      break;

    case 3:
      pm[3] += dnest_randh() * ( range_model[1].inc - range_model[0].inc );
      wrap(&(pm[3]), range_model[0].inc, range_model[1].inc);
      break;

    case 4:
      pm[4] += dnest_randh() * ( range_model[1].opn - range_model[0].opn );
      wrap(&(pm[4]), range_model[0].opn, range_model[1].opn);
      break;

    case 5:
      pm[5] += dnest_randh() * ( range_model[1].A - range_model[0].A );
      wrap(&(pm[5]), range_model[0].A, range_model[1].A);
      break;

    case 6:
      pm[6] += dnest_randh() * ( range_model[1].Ag - range_model[0].Ag );
      wrap(&(pm[6]), range_model[0].Ag, range_model[1].Ag);
      break;

    case 7:
      pm[7] += dnest_randh() * ( range_model[1].k - range_model[0].k );
      wrap(&(pm[7]), range_model[0].k, range_model[1].k);

    case 8:
      pm[8] += dnest_randh() * ( range_model[1].mbh - range_model[0].mbh );
      wrap(&(pm[8]), range_model[0].mbh, range_model[1].mbh);

    case 9:
      pm[9] += dnest_randh() * ( range_model[1].lambda - range_model[0].lambda );
      wrap(&(pm[9]), range_model[0].lambda, range_model[1].lambda);

    case 10:
      pm[10] += dnest_randh() * ( range_model[1].q - range_model[0].q );
      wrap(&(pm[10]), range_model[0].q, range_model[1].q);


    case 11:
      pm[11] += 3.0*dnest_randh();
      wrap(&(pm[11]), -3.0, 0.0);
      break;
    
    case 12:
      pm[12] += 8.0*dnest_randh();
      wrap(&(pm[12]), 2.0, 10.0);
      break;

    case 13:
      pm[13] += 2.0*dnest_randh();
      wrap(&(pm[13]), 0.0, 2.0);
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

void print_particle_line2d(FILE *fp, const void *model)
{
  int i;
  double *pm = (double *)model;

  for(i=0; i<num_params; i++)
  {
    fprintf(fp, "%f ", pm[i] );
  }
  fprintf(fp, "\n");
}

void copy_model_line2d(void *dest, const void *src)
{
  memcpy(dest, src, size_of_modeltype);
}

void *create_model_line2d()
{
  return (void *)malloc(size_of_modeltype);
}

int get_num_params_line2d()
{
  return num_params;
}

void copy_best_model_line2d(const void *bm, const void *bm_std)
{
  memcpy(best_model_line2d, bm, size_of_modeltype);
  memcpy(best_model_std_line2d, bm_std, size_of_modeltype);
}