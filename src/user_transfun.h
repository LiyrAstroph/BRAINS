/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 * \file user_transfun.h
 * \brief header file for user_transfun.c.
 */

#ifndef _BRAINS_USER_TRANSFUN_H

#define _BRAINS_USER_TRANSFUN_H

extern const int num_params_MyTransfun1d;  /*!< number of parameters for 1D model */
extern const int num_params_MyTransfun2d;  /*!< number of parameters for 1D model */

void set_par_range_mytransfun();
void transfun_1d_cal_mytransfun(const void *pm, int flag_save);
void transfun_2d_cal_mytransfun(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save);
void set_par_value_mytransfun_sim(double *pm);
#endif