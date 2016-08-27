/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_interp.h>

#define GRAVITY       6.672e-8
#define SOLAR_MASS    1.989e33
#define C             2.9979e10
#define SEC_PER_YEAR  3.155e7
#define CM_PER_LD    (C*8.64e4)

#define PI            3.141592653589793
#define BRAINS_MAX_STR_LENGTH  (100)

/* MPICH */
extern int thistask, totaltask, namelen;
extern int roottask;
extern char proc_name[MPI_MAX_PROCESSOR_NAME];

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
       tran2d_out_file[BRAINS_MAX_STR_LENGTH];

  char file_dir[BRAINS_MAX_STR_LENGTH];

  int flag_dim;

  int n_con_recon, n_line_recon, n_vel_recon;

  int n_tau;
  double tau_min_set, tau_max_set;

  int n_cloud_per_task, n_vel_per_cloud;

  int flag_save_clouds;

  double InstRes;

  int num_particles;
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
extern double *beta_old_particles;

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

/* BLR model */
typedef struct
{
	double mu;       // in light day;
	double beta;     // 
	double F;        //
	double inc;      
	double opn;	
	double A;
	double Ag;
	double k;
	double mbh;      // in 10e6 solar mass
	double lambda;
	double q;
	//double logse;
}BLRmodel;

extern BLRmodel range_model[2];

/* transfer function / velocity-delay map */
extern double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
extern double dTransTau, dTransV;
extern double rcloud_min_set, rcloud_max_set;

extern double **Fcon_particles;
extern double **Trans1D_particles, **Trans2D_at_veldata_particles; // transfer function 1d for each particle
extern double **clouds_particles;

/* GSL */
extern const gsl_rng_type * gsl_T;
extern gsl_rng * gsl_r;

extern gsl_interp_accel *gsl_acc;
extern gsl_interp  *gsl_linear;

#endif