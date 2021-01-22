/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 * \file user_blr_model.h
 * \brief header file for user_blr_model.c.
 */

#ifndef _BRAINS_USER_BLR_MODEL_H

#define _BRAINS_USER_BLR_MODEL_H

/*
 * define a BLR model struct.
 */
typedef struct
{
  double mu;            /*!< \brief 1. peak radius*/
  double beta;          /*!< \brief 2. fraction of inner radius*/
  double F;             /*!< \brief 3. fraction of outer radius*/
  double inc;           /*!< \brief 4. inclination */
  double opn;           /*!< \brief 5. opening angle */
  double k;             /*!< \brief 6. kappa */
  double gam;           /*!< \brief 7. gamma */
  double xi;            /*!< \brief 8. obscuration */
  double mbh;           /*!< \brief 9. black hole mass */
  double fellip;        /*!< \brief 10.ellipitic orbits */
  double fflow;         /*!< \brief 11.inflow/outflow */
  double sigr_circ;     /*!< \brief 12. */
  double sigthe_circ;   /*!< \brief 13. */
  double sigr_rad;      /*!< \brief 14. */
  double sigthe_rad;    /*!< \brief 15. */
  double theta_rot;     /*!< \brief 16. */
  double sig_turb;      /*!< \brief 17. turbulence velocity */
}MyBLRmodel;

extern const int num_params_MyBLRmodel1d;  /*!< number of parameters for 1D model */
extern const int num_params_MyBLRmodel2d;  /*!< number of parameters for 2D model */

void set_blr_range_mymodel();
void gen_cloud_sample_mymodel(const void *pm, int flag_type, int flag_save);
void set_par_value_mymodel_sim(double *pm);

#ifdef SpecAstro
void set_sa_blr_range_mymodel();
void set_sa_par_value_mymodel_sim(double *pm);
#endif
#endif