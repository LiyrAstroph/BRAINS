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

  int n_con_recon, n_line_recon, n_vel_recon;

  int n_tau;
  double tau_min_set, tau_max_set;

  int n_cloud_per_task, n_vel_per_cloud;

  int flag_save_clouds;

  double InstRes;

  int num_particles;

  char str_par_fix[BRAINS_MAX_STR_LENGTH], str_par_fix_val[BRAINS_MAX_STR_LENGTH];

  int flag_narrowline;
  double flux_narrowline, width_narrowline, shift_narrowline;
}PARSET;
extern PARSET parset;

extern double VelUnit;

extern int n_con_data, n_line_data, n_vel_data;
extern double *Tcon_data, *Fcon_data,  *Fcerrs_data;
extern double *Tline_data, *Fline_data, *Flerrs_data;
extern double *Vline_data, *Fline2d_data, *Flerrs2d_data;
extern double con_scale, line_scale;

extern char dnest_options_file[BRAINS_MAX_STR_LENGTH];

extern int which_parameter_update, which_particle_update;  // which parameter and particle to be updated
extern int which_level_update;
extern int *perturb_accept, *which_parameter_update_prev;
extern double *limits;  // external from dnest

/* continuum reconstruction */
extern double *Tcon, *Fcon, *Fcerrs;
extern double Tcon_min, Tcon_max;
extern double *PSmat, *USmat, *PSmat_data;

/* line reconstruction */
extern double *Fline_at_data;
extern double *Tline, *Fline, *Flerrs;
extern double Tline_min, Tline_max;

/* line reconstruction */
extern double *Fline2d_at_data;
extern double *Fline2d, *Flerrs2d;

/*!
 *  \struct BLRmodel
 *  \brief broad-line region model. 
 */
typedef struct
{
	double mu;       /*!< \brief mean BLR radius, in light day */
	double beta;     /*!< \brief shape parameter */ 
	double F;        /*!< \brief inner edge */
	double inc;      /*!< \brief inclination, in degree, 0-90 */
	double opn;	     /*!< \brief openning angle, in degere, 0-90 */
	double A;        /*!< \brief response coefficient */
	double Ag;       /*!< \brief nonlinear response */
	double k;        /*!< \brief anisotropic emission */
	double mbh;      /*!< \brief black hole mass,  in 10e6 solar mass */
	double lambda;   /*!< \brief orbit parameter */
	double q;        /*!< \brief inflow/outflow */
	double logse;    /*!< \brief systematic errors in continuum and emission line */
}BLRmodel;

extern BLRmodel range_model[2];
extern int *par_fix, npar_fix;
extern double *par_fix_val;

extern int num_params, num_params_blr, num_params_var;

extern double var_range_model[5][2];

/* transfer function / velocity-delay map */
extern double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
extern double dTransTau, dTransV;
extern double rcloud_min_set, rcloud_max_set;

extern double **Fcon_particles, **Fcon_particles_perturb;
extern double *prob_con_particles, *prob_con_particles_perturb;

/* GSL */
extern const gsl_rng_type * gsl_T;
extern gsl_rng * gsl_r;

extern gsl_interp_accel *gsl_acc;
extern gsl_interp  *gsl_linear;

#endif