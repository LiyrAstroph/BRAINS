/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifndef _DNEST_CON_H

#include <stdbool.h>

/*===========================================*/
// users responsible for following struct definitions 

// model

/*==========================================*/
extern int num_params;
extern int size_of_modeltype;

/* best model */
extern void *best_model_thismodel, *best_model_std_thismodel;

/* functions */
void from_prior_thismodel(void *model);
void print_particle_thismodel(FILE *fp, const void *model);
double log_likelihoods_cal_thismodel(const void *model);
double perturb_thismodel(void *model);
void copy_model_thismodel(void *dest, const void *src);
void* create_model_thismodel();
int get_num_params_thismodel();
void copy_best_model_thismodel(const void *bm, const void *bm_std);

void (*print_particle)(FILE *fp, const void *model);
void (*from_prior)(void *model);
double (*log_likelihoods_cal)(const void *model);
double (*perturb)(void *model);
void (*copy_model)(void *dest, const void *src);
void* (*create_model)();
int (*get_num_params)();
void (*copy_best_model)(const void *bm, const void *bm_std);
#endif
