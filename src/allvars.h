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
#include <float.h>
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
#define CM_PER_PC    (3.08567758e+18)
#define CM_PER_MPC   (3.08567758e+24)

#define PI            M_PI
#define BRAINS_MAX_STR_LENGTH  (256)

#define EPS (DBL_MIN)  /* epsilon of the machine as in Matlab */

enum PRIOR_TYPE {GAUSSIAN=1, UNIFORM=2};

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

  int flag_InstRes;
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
  int flag_rng_seed, rng_seed;
  int flag_fixvar;
  int flag_nonlinear;

  int flag_force_update;

  int flag_con_sys_err, flag_line_sys_err;
  
  double temperature;

  int flag_help, flag_end;
  int flag_para_name;
  
  double redshift;
  double linecenter;
#ifdef SA
  char sa_file[BRAINS_MAX_STR_LENGTH];
  int flag_sa_blrmodel, flag_sa_par_mutual;
  double sa_linecenter;
#endif
}PARSET;
extern PARSET parset;

extern double VelUnit, C_Unit;

extern int n_con_data, n_line_data, n_vel_data, n_vel_data_incr, n_vel_data_ext, n_con_max;
extern double *Tcon_data, *Fcon_data,  *Fcerrs_data;
extern double *Tline_data, *Fline_data, *Flerrs_data;
extern double *Vline_data, *Fline2d_data, *Flerrs2d_data;
extern double *Vline_data_ext;
extern double con_scale, line_scale;

extern double con_error_mean, line_error_mean;

extern char dnest_options_file[BRAINS_MAX_STR_LENGTH];

extern int which_parameter_update, which_particle_update;  // which parameter and particle to be updated
extern int which_level_update;
extern double *limits;  // external from dnest

/* continuum reconstruction */
extern int nq;
extern double *Tcon, *Fcon, *Fcerrs;
extern double Tcon_min, Tcon_max;
extern double *PSmat, *PNmat, *USmat, *PSmat_data;
extern double *PNmat_data;
extern double *PCmat_data, *IPCmat_data, *PQmat, *PEmat1, *PEmat2;
extern double *workspace, *workspace_uv;
extern double *var_param, *var_param_std;
extern double *con_q;
extern double logz_con;
extern double *Larr_data, *Larr_rec;
extern double *pow_Tcon_data;

extern double Tspan_data, Tspan_data_con, Tcad_data, Tmed_data;

/* line reconstruction */
extern double *Fline_at_data;
extern double *Tline, *Fline, *Flerrs;
extern double Tline_min, Tline_max;
extern double logz_line;

/* line 2d reconstruction */
extern double *Fline2d_at_data;
extern double *Fline2d, *Flerrs2d;
extern double logz_line2d;

extern int BLRmodel_size;
extern int *par_fix, npar_fix;
extern double *par_fix_val;

extern int num_params_blr_tot;
extern int num_params, num_params_blr, num_params_blr_model, num_params_var, num_params_difftrend, num_params_nlr, num_params_res;
extern int num_params_drw, num_params_trend;
extern int num_params_linecenter;

extern double **blr_range_model, **par_range_model;
extern int *par_prior_model;  /* prior type for model parameters */
extern double **par_prior_gaussian; /* center and std of Gaussian priors */
extern double var_range_model[15][2];
extern double nlr_range_model[3][2];
extern int nlr_prior_model[3];
extern double mass_range[2];
extern double sys_err_line_range[2]; /* range for systematic error of line */
extern double resp_range[2][2];

extern double *instres_epoch, *instres_err_epoch;

/* transfer function / velocity-delay map */
extern double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
extern double rcloud_min_set, rcloud_max_set, time_back_set;

extern double **Fcon_particles, **Fcon_particles_perturb;
extern double *prob_con_particles, *prob_con_particles_perturb;
extern double **con_q_particles, **con_q_particles_perturb;

extern int force_update;
extern double **TransTau_particles, **TransTau_particles_perturb;
extern double **Trans1D_particles, **Trans1D_particles_perturb;
extern double **Trans2D_at_veldata_particles, **Trans2D_at_veldata_particles_perturb;
extern double **Fline_at_data_particles, **Fline_at_data_particles_perturb;

extern double *clouds_tau, *clouds_weight, *clouds_vel;

extern FILE *fcloud_out;
extern int icr_cloud_save;

extern int postprc;

/* GSL */
extern const gsl_rng_type * gsl_T;
extern gsl_rng * gsl_r;

extern gsl_interp_accel *gsl_acc;
extern gsl_interp  *gsl_linear;


#ifdef SA
extern double PhaseFactor;

extern int num_params_rm;
extern int num_params_sa, num_params_sa_blr_model, num_params_sa_extpar, num_params_sa_blr;
extern int n_epoch_sa_data, n_vel_sa_data, n_base_sa_data;
extern double *vel_sa_data, *base_sa_data, *Fline_sa_data, *Flerrs_sa_data, *phase_sa_data, *pherrs_sa_data;
extern double *wave_sa_data;

extern double sa_flux_norm;

extern double **phase_sa_particles, **Fline_sa_particles;
extern double **phase_sa_particles_perturb, **Fline_sa_particles_perturb;
extern double *phase_sa, *Fline_sa;

extern double *clouds_alpha, *clouds_beta;

extern double **sa_extpar_range, **sa_blr_range_model;

extern int SABLRmodel_size;

extern int *idx_sa_par_mutual, *idx_rm_par_mutual;

extern double *prob_sa_particles, *prob_sa_particles_perturb;

extern double *workspace_phase;
#endif
#endif