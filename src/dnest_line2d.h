/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file dnest_line2d.h
 *  \brief header file.
 */

#ifndef _DNEST_LINE2D_H
#define _DNEST_LINE2D_H

#include <stdbool.h>

/* functions */
void from_prior_line2d(void *model);
void print_particle_line2d(FILE *fp, const void *model);
void read_particle_line2d(FILE *fp, void *model);
void restart_action_2d(int iflag);
void accept_action_2d();
void kill_action_2d(int i, int i_copy);
double log_likelihoods_cal_line2d(const void *model);
double log_likelihoods_cal_initial_line2d(const void *model);
double log_likelihoods_cal_restart_line2d(const void *model);
double perturb_line2d(void *model);
double log_likelihoods_cal_line2d_exam(const void *model);

#endif