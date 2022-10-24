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

extern char **BLRmodel_name;
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
  double mbh;      /*!< \brief 7.  black hole mass,  in 10e6 solar mass */
  double lambda;   /*!< \brief 8.  orbit parameter */
  double q;        /*!< \brief 9.  inflow/outflow */
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
  double mbh;      /*!< \brief 7.  black hole mass,  in 10e6 solar mass */
  double sigr;     /*!< \brief 8.  orbit parameter */
  double sigtheta; /*!< \brief 9.  inflow/outflow */
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
  double mbh;      /*!< \brief 7.  black hole mass,  in 10e6 solar mass */
  double xi;       /*!< \brief 8.  orbit parameter */
  double q;        /*!< \brief 9.  inflow/outflow */
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
  double mbh;           /*!< \brief 10.black hole mass */
  double fellip;        /*!< \brief 11.elliptic orbits */
  double fflow;         /*!< \brief 12.inflow/outflow */
  double sigr_circ;     /*!< \brief 13. */
  double sigthe_circ;   /*!< \brief 14. */
  double sigr_rad;      /*!< \brief 15. */
  double sigthe_rad;    /*!< \brief 16. */
  double theta_rot;     /*!< \brief 17. */
  double sig_turb;      /*!< \brief 18.turbulence velocity */
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
  double mbh;           /*!< \brief 9. black hole mass */
  double fellip;        /*!< \brief 10.elliptic orbits */
  double fflow;         /*!< \brief 11.inflow/outflow */
  double sigr_circ;     /*!< \brief 12. */
  double sigthe_circ;   /*!< \brief 13. */
  double sigr_rad;      /*!< \brief 14. */
  double sigthe_rad;    /*!< \brief 15. */
  double theta_rot;     /*!< \brief 16. */
  double sig_turb;      /*!< \brief 17. turbulence velocity */
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

  double mbh;           /*!< \brief 14.black hole mass */
  double fellip;        /*!< \brief 15.elliptic orbits */
  double fflow;         /*!< \brief 16.inflow/outflow */
  double sigr_circ;     /*!< \brief 17. */
  double sigthe_circ;   /*!< \brief 18. */
  double sigr_rad;      /*!< \brief 19. */
  double sigthe_rad;    /*!< \brief 20. */
  double theta_rot;     /*!< \brief 21. */
  double fellip_un;     /*!< \brief 22.elliptic orbits */
  double fflow_un;      /*!< \brief 23.inflow/outflow */
  double sig_turb;      /*!< \brief 24. */
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
  double mbh;           /*!< \brief 13 */
}BLRmodel8;

typedef struct
{
  double mu;       /*!< \brief 1.  mean BLR radius, in light day */
  double beta;     /*!< \brief 2.  shape parameter */ 
  double F;        /*!< \brief 3.  inner edge */
  double inc;      /*!< \brief 4.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 5.  openning angle, in degere, 0-90 */
  double mbh;      /*!< \brief 6.  black hole mass,  in 10e6 solar mass */ 
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
}SAExtPar;

#endif

#endif