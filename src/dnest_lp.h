/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_lp.h
 *  \brief header file.
 */

#ifndef _DNEST_LP_H
#define _DNEST_LP_H

#include <stdbool.h>

/* functions */
void from_prior_lp(void *model);
void print_particle_lp(FILE *fp, const void *model);
void read_particle_lp(FILE *fp, void *model);
void restart_action_lp(int iflag);
void accept_action_lp();
void kill_action_lp(int i, int i_copy);
double log_likelihoods_cal_lp(const void *model);
double log_likelihoods_cal_initial_lp(const void *model);
double log_likelihoods_cal_restart_lp(const void *model);
double perturb_lp(void *model);
double log_likelihoods_cal_lp_exam(const void *model);

#endif