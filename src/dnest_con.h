/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
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
double log_likelihoods_cal_con(const void *model);
double log_likelihoods_cal_initial_con(const void *model);
double log_likelihoods_cal_restart_con(const void *model);
double perturb_con(void *model);
void restart_action_con(int iflag);

double log_likelihoods_cal_con_exam(const void *model);

#endif
