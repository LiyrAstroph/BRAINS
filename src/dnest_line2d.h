/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
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

/* size of model type, defined in dnest */
extern int size_of_modeltype;

/* functions */
void from_prior_line2d(void *model);
void print_particle_line2d(FILE *fp, const void *model);
double log_likelihoods_cal_line2d(const void *model);
double log_likelihoods_cal_initial_line2d(const void *model);
double log_likelihoods_cal_restart_line2d(const void *model);
//double perturb_line2d(void *model);
int get_num_params_line2d();

double perturb_line2d_model1(void *model);
double perturb_line2d_model3(void *model);
double perturb_line2d_model5(void *model);
double perturb_line2d_model6(void *model);

double log_likelihoods_cal_line2d_exam(const void *model);


void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*log_likelihoods_cal_initial)(const void *model);
double (*log_likelihoods_cal_restart)(const void *model);
double (*perturb)(void *model);
int (*get_num_params)();
void (*restart_clouds)(int iflag);
#endif