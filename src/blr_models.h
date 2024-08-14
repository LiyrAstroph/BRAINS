/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file blr_models.h
 *  \brief header file of blr_models.c.
 */

#ifndef _BLR_MODELS_H
#define _BLR_MODELS_H

double (*theta_sample)(double gam, double Lopn_cos1, double Lopn_cos2);
double theta_sample_outer(double gam, double Lopn_cos1, double Lopn_cos2);
double theta_sample_inner(double gam, double Lopn_cos1, double Lopn_cos2);
double eta_func_powerlaw(double eta0, double eta1, double alpha, double r);
double eta_func_step(double eta0, double eta1, double alpha, double r);
double (*eta_func)(double eta0, double eta1, double alpha, double r);

extern char **BLRmodel_name;
extern char **BLRmodel_sa_name;
extern char *BLRmodel1_name[];
extern char *BLRmodel2_name[];
extern char *BLRmodel3_name[];
extern char *BLRmodel4_name[];
extern char *BLRmodel5_name[];
extern char *BLRmodel6_name[];
extern char *BLRmodel7_name[];
extern char *BLRmodel8_name[];
extern char *BLRmodel9_name[];

/*!
 *  \struct BLRmodel1
 *  \brief broad-line region model 1. 
 */
typedef struct
{
  double mu;       /*!< \brief 1.  mean BLR radius, in light day */
  double beta;     /*!< \brief 2.  shape parameter */ 
  double F;        /*!< \brief 3.  inner edge */
  double inc;      /*!< \brief 4.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 5.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 6.  anisotropic emission */
  double eta0;     /*!< \brief 7. */
  double eta1;     /*!< \brief 8. */
  double eta_alpha;/*!< \brief 9. */
  double mbh;      /*!< \brief 10.  black hole mass,  in 10e6 solar mass */
  double lambda;   /*!< \brief 11.  orbit parameter */
  double q;        /*!< \brief 12.  inflow/outflow */
}BLRmodel1;

/*!
 *  \struct BLRmodel2
 *  \brief broad-line region model 2. 
 */
typedef struct
{
  double mu;       /*!< \brief 1.  mean BLR radius, in light day */
  double beta;     /*!< \brief 2.  shape parameter */ 
  double F;        /*!< \brief 3.  inner edge */
  double inc;      /*!< \brief 4.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 5.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 6.  anisotropic emission */
  double eta0;     /*!< \brief 7. */
  double eta1;     /*!< \brief 8. */
  double eta_alpha;/*!< \brief 9. */
  double mbh;      /*!< \brief 10.  black hole mass,  in 10e6 solar mass */
  double sigr;     /*!< \brief 11.  orbit parameter */
  double sigtheta; /*!< \brief 12.  inflow/outflow */
}BLRmodel2;

/*!
 *  \struct BLRmodel3
 *  \brief broad-line region model 3. 
 */
typedef struct
{
  double Rin;      /*!< \brief 1.  mean BLR radius, in light day */
  double F;        /*!< \brief 2.  inner edge */
  double alpha;    /*!< \brief 3.  shape parameter */ 
  double inc;      /*!< \brief 4.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 5.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 6.  anisotropic emission */
  double eta0;     /*!< \brief 7. */
  double eta1;     /*!< \brief 8. */
  double eta_alpha;/*!< \brief 9. */
  double mbh;      /*!< \brief 10.  black hole mass,  in 10e6 solar mass */
  double xi;       /*!< \brief 11.  orbit parameter */
  double q;        /*!< \brief 12.  inflow/outflow */
}BLRmodel3;

/*!
 *  \struct BLRmodel4
 *  \brief broad-line region model 4. 
 */
typedef BLRmodel3 BLRmodel4;

/*!
 *  \struct BLRmodel5
 *  \brief broad-line region model 5. 
 */
typedef struct 
{
  double mu;            /*!< \brief 1. peak radius*/
  double Fin;           /*!< \brief 2. fraction of inner radius*/
  double Fout;          /*!< \brief 3. fraction of outer radius*/
  double alpha;         /*!< \brief 4. slope*/
  double inc;           /*!< \brief 5. inclination */
  double opn;           /*!< \brief 6. opening angle */
  double k;             /*!< \brief 7. kappa */
  double gam;           /*!< \brief 8. gamma */
  double xi;            /*!< \brief 9. obscuration */
  double eta0;          /*!< \brief 10. */
  double eta1;          /*!< \brief 11. */
  double eta_alpha;     /*!< \brief 12. */
  double mbh;           /*!< \brief 13.black hole mass */
  double fellip;        /*!< \brief 14.elliptic orbits */
  double fflow;         /*!< \brief 15.inflow/outflow */
  double sigr_circ;     /*!< \brief 16. */
  double sigthe_circ;   /*!< \brief 17. */
  double sigr_rad;      /*!< \brief 18. */
  double sigthe_rad;    /*!< \brief 19. */
  double theta_rot;     /*!< \brief 20. */
  double sig_turb;      /*!< \brief 21.turbulence velocity */
}BLRmodel5;

/*!
 *  \struct BLRmodel6
 *  \brief broad-line region model 6, Pancoast 2014's model. 
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
  double eta0;          /*!< \brief 9. */
  double eta1;          /*!< \brief 10. */
  double eta_alpha;     /*!< \brief 11. */
  double mbh;           /*!< \brief 12. black hole mass */
  double fellip;        /*!< \brief 13.elliptic orbits */
  double fflow;         /*!< \brief 14.inflow/outflow */
  double sigr_circ;     /*!< \brief 15. */
  double sigthe_circ;   /*!< \brief 16. */
  double sigr_rad;      /*!< \brief 17. */
  double sigthe_rad;    /*!< \brief 18. */
  double theta_rot;     /*!< \brief 19. */
  double sig_turb;      /*!< \brief 20. turbulence velocity */
}BLRmodel6;


/*!
 *  \struct BLRmodel7
 *  \brief broad-line region model 7, shadowed model in Li et al. 2018. 
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

  double fsh;           /*!< \brief 9. fraction of two regions */
  double mu_un;         /*!< \brief 10.peak radius*/
  double beta_un;       /*!< \brief 11.fraction of inner radius*/
  double F_un;          /*!< \brief 12.fraction of outer radius*/
  double opn_un;        /*!< \brief 13.opening angle */

  double eta0;          /*!< \brief 14. */
  double eta1;          /*!< \brief 15. */
  double eta_alpha;     /*!< \brief 16. */

  double mbh;           /*!< \brief 17.black hole mass */
  double fellip;        /*!< \brief 18.elliptic orbits */
  double fflow;         /*!< \brief 19.inflow/outflow */
  double sigr_circ;     /*!< \brief 20. */
  double sigthe_circ;   /*!< \brief 21. */
  double sigr_rad;      /*!< \brief 22. */
  double sigthe_rad;    /*!< \brief 23. */
  double theta_rot;     /*!< \brief 24. */
  double fellip_un;     /*!< \brief 25.elliptic orbits */
  double fflow_un;      /*!< \brief 26.inflow/outflow */
  double sig_turb;      /*!< \brief 27. */
}BLRmodel7;


/*!
 *  \struct BLRmodel8
 *  \brief broad-line region model 8, disk wind model
 */
typedef struct
{
  double theta_min;     /*!< \brief 1 */
  double dtheta_max;    /*!< \brief 2 */
  double r_min;         /*!< \brief 3 */
  double fr_max;        /*!< \brief 4 */
  double gamma;         /*!< \brief 5 */
  double alpha;         /*!< \brief 6 */
  double lamda;         /*!< \brief 7 */
  double k;         /*!< \brief 8 */
  double xi;            /*!< \brief 9 */
  double Rv;            /*!< \brief 10 */ 
  double Rblr;          /*!< \brief 11 */
  double inc;           /*!< \brief 12 */
  double eta0;          /*!< \brief 13. */
  double eta1;          /*!< \brief 14. */
  double eta_alpha;     /*!< \brief 15. */
  double mbh;           /*!< \brief 16 */
}BLRmodel8;

typedef struct
{
  double mu;       /*!< \brief 1.  mean BLR radius, in light day */
  double beta;     /*!< \brief 2.  shape parameter */ 
  double F;        /*!< \brief 3.  inner edge */
  double inc;      /*!< \brief 4.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 5.  openning angle, in degere, 0-90 */
  double eta0;     /*!< \brief 6. */
  double eta1;     /*!< \brief 7. */
  double eta_alpha;/*!< \brief 8. */
  double mbh;      /*!< \brief 9.  black hole mass,  in 10e6 solar mass */ 
}BLRmodel9;

#ifdef SpecAstro

typedef BLRmodel1 SABLRmodel1;
typedef BLRmodel2 SABLRmodel2;
typedef BLRmodel3 SABLRmodel3;
typedef BLRmodel4 SABLRmodel4;
typedef BLRmodel5 SABLRmodel5;
typedef BLRmodel6 SABLRmodel6;
typedef BLRmodel7 SABLRmodel7;
typedef BLRmodel8 SABLRmodel8;
typedef BLRmodel9 SABLRmodel9;

typedef struct 
{
  double DA; /* angular size distance */
  double PA; /* position angle */
  double FA; /* line flux scaling factor */
  double CO; /* line center offset */
  double AlphaC; /* offset wrt to continuum */
  double BetaC;  /* offset wrt to continuum */
}SAExtPar;

#endif

#endif