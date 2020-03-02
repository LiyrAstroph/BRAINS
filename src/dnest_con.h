/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_con.h
 *  \brief header file.
 */

#ifndef _DNEST_CON_H

#define _DNEST_CON_H

#include <stdbool.h>

/* functions */
void from_prior_con(void *model);
void print_particle_con(FILE *fp, const void *model);
void read_particle_con(FILE *fp, void *model);
void restart_action_con(int iflag);
void accept_action_con();
void kill_action_con(int i, int i_copy);
double log_likelihoods_cal_con(const void *model);
double log_likelihoods_cal_initial_con(const void *model);
double log_likelihoods_cal_restart_con(const void *model);
double perturb_con(void *model);
double log_likelihoods_cal_con_exam(const void *model);

#endif
