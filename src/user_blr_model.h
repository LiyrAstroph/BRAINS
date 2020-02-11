/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 * \file user_blr_model.h
 * \brief header file for user_blr_model.c.
 */

#ifndef _BRAINS_USER_BLR_MODEL_H

#define _BRAINS_USER_BLR_MODEL_H

extern const int num_params_MyBLRmodel1d;  /*!< number of parameters for 1D model */
extern const int num_params_MyBLRmodel2d;  /*!< number of parameters for 2D model */

void set_blr_range_mymodel();
void gen_cloud_sample_mymodel(const void *pm, int flag_type, int flag_save);
void set_par_value_mymodel_sim(double *pm);
#endif