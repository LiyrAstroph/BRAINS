/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifdef SpecAstro

#include "brains.h"

#ifndef _DNEST_SARM_H
#define _DNEST_SARM_H

#include <stdbool.h>

/* functions */
void from_prior_sarm(void *model);
void print_particle_sarm(FILE *fp, const void *model);
void read_particle_sarm(FILE *fp, void *model);
void restart_action_sarm(int iflag);
void accept_action_sarm();
void kill_action_sarm(int i, int i_copy);
double log_likelihoods_cal_sarm(const void *model);
double log_likelihoods_cal_initial_sarm(const void *model);
double log_likelihoods_cal_restart_sarm(const void *model);
double perturb_sarm(void *model);
double log_likelihoods_cal_sarm_exam(const void *model);

#endif 

#endif