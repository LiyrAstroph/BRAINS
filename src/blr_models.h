/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifndef _BLR_MODELS_H
#define _BLR_MODELS_H

/*!
 *  \struct BLRmodel1
 *  \brief broad-line region model. 
 */
typedef struct
{
  double A;        /*!< \brief 1.  response coefficient */
  double Ag;       /*!< \brief 2.  nonlinear response */
  double mu;       /*!< \brief 3.  mean BLR radius, in light day */
  double beta;     /*!< \brief 4.  shape parameter */ 
  double F;        /*!< \brief 5.  inner edge */
  double inc;      /*!< \brief 6.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 7.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 8.  anisotropic emission */
  double mbh;      /*!< \brief 9.  black hole mass,  in 10e6 solar mass */
  double lambda;   /*!< \brief 10. orbit parameter */
  double q;        /*!< \brief 11. inflow/outflow */
  double logse;    /*!< \brief 12. systematic errors in continuum and emission line */
}BLRmodel1;

typedef struct
{
  double A;        /*!< \brief 1.  response coefficient */
  double Ag;       /*!< \brief 2.  nonlinear response */
  double mu;       /*!< \brief 3.  mean BLR radius, in light day */
  double beta;     /*!< \brief 4.  shape parameter */ 
  double F;        /*!< \brief 5.  inner edge */
  double inc;      /*!< \brief 6.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 7.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 8.  anisotropic emission */
  double mbh;      /*!< \brief 9.  black hole mass,  in 10e6 solar mass */
  double sigr;     /*!< \brief 10. orbit parameter */
  double sigtheta; /*!< \brief 11. inflow/outflow */
  double logse;    /*!< \brief 12. systematic errors in continuum and emission line */
}BLRmodel2;

typedef struct
{
  double A;        /*!< \brief 1.  response coefficient */
  double Ag;       /*!< \brief 2.  nonlinear response */
  double alpha;    /*!< \brief 3.  shape parameter */ 
  double Rin;      /*!< \brief 4.  mean BLR radius, in light day */
  double F;        /*!< \brief 5.  inner edge */
  double inc;      /*!< \brief 6.  inclination, in degree, 0-90 */
  double opn;      /*!< \brief 7.  openning angle, in degere, 0-90 */
  double k;        /*!< \brief 8.  anisotropic emission */
  double mbh;      /*!< \brief 9.  black hole mass,  in 10e6 solar mass */
  double xi;       /*!< \brief 10. orbit parameter */
  double q;        /*!< \brief 11. inflow/outflow */
  double logse;    /*!< \brief 12. systematic errors in continuum and emission line */
}BLRmodel3;


typedef BLRmodel3 BLRmodel4;

typedef struct 
{
  double A;             /*!< \brief 1. line response*/    
  double Ag;            /*!< \brief 2. no-linearity */
  double mu;            /*!< \brief 3. peak radius*/
  double Fin;           /*!< \brief 4. fraction of inner radius*/
  double Fout;          /*!< \brief 5. fraction of outer radius*/
  double alpha;         /*!< \brief 6. slope*/
  double inc;           /*!< \brief 7. inclination */
  double opn;           /*!< \brief 8. opening angle */
  double k;             /*!< \brief 9. kappa */
  double gam;           /*!< \brief 10.gamma */
  double xi;            /*!< \brief 11.obscuration */
  double mbh;           /*!< \brief 12.black hole mass */
  double fellip;        /*!< \brief 13.ellipitic orbits */
  double fflow;         /*!< \brief 14.inflow/outflow */
  double sigr_circ;     /*!< \brief 15. */
  double sigthe_circ;   /*!< \brief 16. */
  double sigr_rad;      /*!< \brief 17. */
  double sigthe_rad;    /*!< \brief 18. */
  double theta_rot;     /*!< \brief 19. */
  double sig_turb;      /*!< \brief 20. turbulence velocity */
  double logse;         /*!< \brief 21. */
}BLRmodel5;

/* Pancoast 2014's model */
typedef struct 
{
  double A;             /*!< \brief 1. line response*/    
  double Ag;            /*!< \brief 2. no-linearity */
  double mu;            /*!< \brief 3. peak radius*/
  double beta;          /*!< \brief 4. fraction of inner radius*/
  double F;             /*!< \brief 5. fraction of outer radius*/
  double inc;           /*!< \brief 6. inclination */
  double opn;           /*!< \brief 7. opening angle */
  double k;             /*!< \brief 8. kappa */
  double gam;           /*!< \brief 9. gamma */
  double xi;            /*!< \brief 10.obscuration */
  double mbh;           /*!< \brief 11.black hole mass */
  double fellip;        /*!< \brief 12.ellipitic orbits */
  double fflow;         /*!< \brief 13.inflow/outflow */
  double sigr_circ;     /*!< \brief 14. */
  double sigthe_circ;   /*!< \brief 15. */
  double sigr_rad;      /*!< \brief 16. */
  double sigthe_rad;    /*!< \brief 17. */
  double theta_rot;     /*!< \brief 18. */
  double sig_turb;      /*!< \brief 19. turbulence velocity */
  double logse;         /*!< \brief 20. */
}BLRmodel6;


/* shadowed model */
typedef struct 
{
  double A;             /*!< \brief 1. line response*/    
  double Ag;            /*!< \brief 2. no-linearity */
  double mu;            /*!< \brief 3. peak radius*/
  double beta;          /*!< \brief 4. fraction of inner radius*/
  double F;             /*!< \brief 5. fraction of outer radius*/
  double inc;           /*!< \brief 6. inclination */
  double opn;           /*!< \brief 7. opening angle */
  double k;             /*!< \brief 8. kappa */
  double gam;           /*!< \brief 9. gamma */
  double xi;            /*!< \brief 10.obscuration */

  double fsh;           /*!< 11 */
  double mu_un;         /*!< \brief 12 peak radius*/
  double beta_un;       /*!< \brief 13 fraction of inner radius*/
  double F_un;          /*!< \brief 14 fraction of outer radius*/
  double opn_un;        /*!< \brief 15 opening angle */

  double mbh;           /*!< \brief 16.black hole mass */
  double fellip;        /*!< \brief 17.ellipitic orbits */
  double fflow;         /*!< \brief 18.inflow/outflow */
  double sigr_circ;     /*!< \brief 19. */
  double sigthe_circ;   /*!< \brief 20. */
  double sigr_rad;      /*!< \brief 21. */
  double sigthe_rad;    /*!< \brief 22. */
  double theta_rot;     /*!< \brief 23. */
  double fellip_un;        /*!< \brief 24.ellipitic orbits */
  double fflow_un;         /*!< \brief 25.inflow/outflow */
  double sig_turb;      /*!< \brief 26. */
  double logse;         /*!< \brief 27. */
}BLRmodel7;

#endif