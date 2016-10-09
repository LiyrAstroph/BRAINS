/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifndef _DNEST_LINE2D_H
#define _DNEST_LINE2D_H

#include <stdbool.h>


/*==========================================*/
//extern int num_params, num_params_blr, num_params_var;
extern int size_of_modeltype;

/* best model */
extern void *best_model_line2d, *best_model_std_line2d;

/* functions */
void from_prior_line2d(void *model);
void print_particle_line2d(FILE *fp, const void *model);
double log_likelihoods_cal_line2d(const void *model);
double perturb_line2d(void *model);
void copy_model_line2d(void *dest, const void *src);
void* create_model_line2d();
int get_num_params_line2d();
void copy_best_model_line2d(const void *bm, const void *bm_std);

void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*perturb)(void *model);
void (*copy_model)(void *dest, const void *src);
void* (*create_model)();
int (*get_num_params)();
void (*copy_best_model)(const void *bm, const void *bm_std);

#endif