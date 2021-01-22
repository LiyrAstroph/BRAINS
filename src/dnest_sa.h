/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_sa.h
 *  \brief header file.
 */

#ifdef SpecAstro

#ifndef _DNEST_SA_H
#define _DNEST_SA_H

#include <stdbool.h>

/* functions */
void from_prior_sa(void *model);
void print_particle_sa(FILE *fp, const void *model);
void read_particle_sa(FILE *fp, void *model);
void restart_action_sa(int iflag);
void accept_action_sa();
void kill_action_sa(int i, int i_copy);
double log_likelihoods_cal_sa(const void *model);
double log_likelihoods_cal_initial_sa(const void *model);
double log_likelihoods_cal_restart_sa(const void *model);
double perturb_sa(void *model);
double log_likelihoods_cal_sa_exam(const void *model);

#endif

#endif