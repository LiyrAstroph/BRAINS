/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_sa.h
 *  \brief header file.
 */

#ifdef SA

#ifndef _DNEST_SA_H
#define _DNEST_SA_H

#include <stdbool.h>

/* functions */
void from_prior_sa2d(void *model);
void print_particle_sa2d(FILE *fp, const void *model);
void read_particle_sa2d(FILE *fp, void *model);
void restart_action_sa2d(int iflag);
void accept_action_sa2d();
void kill_action_sa2d(int i, int i_copy);
double log_likelihoods_cal_sa2d(const void *model);
double log_likelihoods_cal_initial_sa2d(const void *model);
double log_likelihoods_cal_restart_sa2d(const void *model);
double perturb_sa2d(void *model);
double log_likelihoods_cal_sa2d_exam(const void *model);

#endif

#endif