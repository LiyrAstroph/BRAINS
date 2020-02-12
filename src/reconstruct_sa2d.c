/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_sa2d.c
 *  \brief reconstruct sa and 2d RM data.
 */

#ifdef SA

#include "brains.h"

void *best_model_sa2d;      /*!< best model */
void *best_model_std_sa2d;  /*!< standard deviation of the best model */

void reconstruct_sa2d()
{
  int i, argc=0;
  char **argv;
}

void reconstruct_sa2d_init()
{
  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSSA2D");

  /* cloud sample related */
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));
}

#endif