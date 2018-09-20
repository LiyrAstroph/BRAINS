/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file allvars.h
 *  \brief header file of allvars.c.
 */

#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_interp.h>

/*! \def GRAVITY
 *  \brief Gravitational constant. 

 *  \def SOLAR_MASS
 *  \brief Solar mass.

 *  \def C
 *  \brief speed of light.

 *  \def SEC_PER_YEAR
 *  \brief seconds per year.

 *  \def CM_PER_LD
 *  \brief light day in centimeter.

 *  \def PI
 *  \brief Pi.

 *  \def BRAINS_MAX_STR_LENGTH
 *  \brief maximum string length.
 */
#define GRAVITY       6.672e-8
#define SOLAR_MASS    1.989e33
#define C             2.9979e10
#define SEC_PER_YEAR  3.155e7
#define CM_PER_LD    (C*8.64e4)

#define PI            M_PI
#define BRAINS_MAX_STR_LENGTH  (100)

/* variables for MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

/*! \struct PARSET
 *  \brief the configuration parameters.
 */
typedef struct 
{
  char param_file[BRAINS_MAX_STR_LENGTH];
  char continuum_file[BRAINS_MAX_STR_LENGTH], 
       line_file[BRAINS_MAX_STR_LENGTH], 
       line2d_file[BRAINS_MAX_STR_LENGTH],
       pcon_out_file[BRAINS_MAX_STR_LENGTH],
       pline_out_file[BRAINS_MAX_STR_LENGTH],
       pline2d_out_file[BRAINS_MAX_STR_LENGTH],
       pline2d_data_out_file[BRAINS_MAX_STR_LENGTH],
       cloud_out_file[BRAINS_MAX_STR_LENGTH],
       tran_out_file[BRAINS_MAX_STR_LENGTH],
       tran2d_out_file[BRAINS_MAX_STR_LENGTH],
       tran2d_data_out_file[BRAINS_MAX_STR_LENGTH];

  char file_dir[BRAINS_MAX_STR_LENGTH];

  int flag_dim;

  int flag_trend, flag_trend_diff;
  
  int n_con_recon, n_line_recon, n_vel_recon;

  int n_tau;
  double rcloud_max, time_back;

  int flag_blrmodel;

  int n_cloud_per_task, n_vel_per_cloud;

  int flag_save_clouds;

  double InstRes, InstRes_err;
  char file_instres[BRAINS_MAX_STR_LENGTH];

  int num_particles;

  char str_par_fix[BRAINS_MAX_STR_LENGTH], str_par_fix_val[BRAINS_MAX_STR_LENGTH];

  int flag_narrowline;
  double flux_narrowline, width_narrowline, shift_narrowline;
  double flux_narrowline_err, width_narrowline_err, shift_narrowline_err;
  
  int flag_linecenter;
  double linecenter_err;

  int flag_postprc;
  int flag_restart;
  int flag_sample_info;
  int flag_temp;
  int flag_exam_prior;

  int flag_fixvar;

  int flag_force_update;
  
  double temperature;
}PARSET;
extern PARSET parset;

extern double VelUnit, C_Unit;

extern int n_con_data, n_line_data, n_vel_data;
extern double *Tcon_data, *Fcon_data,  *Fcerrs_data;
extern double *Tline_data, *Fline_data, *Flerrs_data;
extern double *Vline_data, *Fline2d_data, *Flerrs2d_data;
extern double con_scale, line_scale;

extern double con_error_mean, line_error_mean;

extern char dnest_options_file[BRAINS_MAX_STR_LENGTH];

extern int which_parameter_update, which_particle_update;  // which parameter and particle to be updated
extern int which_level_update;
extern unsigned long long int which_mcmc_steps;//mcmc steps 
extern int *perturb_accept, *which_parameter_update_prev;
extern double *limits;  // external from dnest

/* continuum reconstruction */
extern int nq;
extern double *Tcon, *Fcon, *Fcerrs;
extern double Tcon_min, Tcon_max;
extern double *PSmat, *PNmat, *USmat, *PSmat_data;
extern double *PNmat_data;
extern double *PCmat_data, *IPCmat_data, *PQmat, *PEmat1, *PEmat2;
extern double *workspace;
extern double *var_param, *var_param_std;
extern double *con_q;

/* line reconstruction */
extern double *Fline_at_data;
extern double *Tline, *Fline, *Flerrs;
extern double Tline_min, Tline_max;

/* line reconstruction */
extern double *Fline2d_at_data;
extern double *Fline2d, *Flerrs2d;

extern int num_params_radial_samp;
extern int *params_radial_samp;

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
	double opn;	     /*!< \brief 7.  openning angle, in degere, 0-90 */
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

extern int BLRmodel_size;
extern int *par_fix, npar_fix;
extern double *par_fix_val;

extern int num_params, num_params_blr, num_params_blr_model, num_params_var, num_params_difftrend, num_params_nlr, num_params_res;
extern int num_params_linecenter;

extern double **blr_range_model, **par_range_model;
extern double var_range_model[7][2];
extern double nlr_range_model[3][2];

extern double *instres_epoch, *instres_err_epoch;

/* transfer function / velocity-delay map */
extern double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
extern double dTransTau, dTransV;
extern double rcloud_min_set, rcloud_max_set;

extern double **Fcon_particles, **Fcon_particles_perturb;
extern double *prob_con_particles, *prob_con_particles_perturb;
extern double **con_q_particles, **con_q_particles_perturb;

extern int force_update;
extern double **clouds_particles_perturb, **clouds_particles;
extern double **Trans1D_particles, **Trans1D_particles_perturb;
extern double **Trans2D_at_veldata_particles, **Trans2D_at_veldata_particles_perturb;
extern double **Fline_at_data_particles, **Fline_at_data_particles_perturb;
extern double *prob_line_particles, *prob_line_particles_perturb;

extern double tau_max, tau_min;
extern double *tmp_tau, *tmp_weight, *tmp_vel;

extern int postprc;

/* GSL */
extern const gsl_rng_type * gsl_T;
extern gsl_rng * gsl_r;

extern gsl_interp_accel *gsl_acc;
extern gsl_interp  *gsl_linear;

#endif