/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file allvars.c
 *  \brief define all variables.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "allvars.h"

/* MPICH */
int thistask, totaltask, namelen;
int roottask = 0;
char proc_name[MPI_MAX_PROCESSOR_NAME];

PARSET parset;


double VelUnit, C_Unit;    /*!< Velocity unit = \f$\sqrt{\frac{G\cdot 10^6M_\odot}{1~\rm light-day}}/10^5\f$km/s, 
                        with \f$ M_\bullet = \f$ and \f$r=1\f$light day.*/

/* Data */
int n_con_data;    /*!< number of continuum data points */
int n_line_data;   /*!< number of emission data points along time axis */
int n_vel_data;    /*!< number of emission data points along velocity axis */
double *Tcon_data, *Fcon_data,  *Fcerrs_data;
double *Tline_data, *Fline_data, *Flerrs_data;
double *Vline_data, *Fline2d_data, *Flerrs2d_data;
double con_scale, line_scale;

double con_error_mean, line_error_mean;

char dnest_options_file[BRAINS_MAX_STR_LENGTH];

int which_parameter_update;  /*!< which parameter being updated */
int which_particle_update;   /*!< which particle being updated */
int which_level_update;      /*!< which level of the particle */

/* continuum reconstruction */
int nq;
double *Tcon, *Fcon, *Fcerrs;
double Tcon_min, Tcon_max;
double *PSmat, *PNmat, *USmat, *PSmat_data;
double *PNmat_data;
double *PCmat_data, *IPCmat_data, *PQmat, *PEmat1, *PEmat2;
double *workspace;
double *var_param, *var_param_std;
double *con_q;
double logz_con;


/* line reconstruction */
double *Fline_at_data;
double *Tline, *Fline, *Flerrs;
double Tline_min, Tline_max;
double logz_line;

/* line 2d reconstruction */
double *Fline2d_at_data;
double *Fline2d, *Flerrs2d;
double logz_line2d;

int num_params_radial_samp;
int *params_radial_samp;

// BLR
int BLRmodel_size;
int *par_fix, npar_fix;
double *par_fix_val;
double **blr_range_model, **par_range_model;
double nlr_range_model[3][2];
double mass_range[2];

int num_params, num_params_blr, num_params_blr_model, num_params_var, num_params_difftrend, num_params_nlr, num_params_res;
int num_params_linecenter;

// continuum variation
double var_range_model[7][2]; /*!< define the range of variability parameters */

// epoch dependent spectral broadening
double *instres_epoch, *instres_err_epoch;

/* transfer function / velocity-delay map */
double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
double dTransTau, dTransV;
double rcloud_min_set, rcloud_max_set;

double **Fcon_particles, **Fcon_particles_perturb;
double *prob_con_particles, *prob_con_particles_perturb;
double **con_q_particles, **con_q_particles_perturb;

int force_update;
double **clouds_particles_perturb, **clouds_particles;
double **Trans1D_particles, **Trans1D_particles_perturb;
double **Trans2D_at_veldata_particles, **Trans2D_at_veldata_particles_perturb;
double **Fline_at_data_particles, **Fline_at_data_particles_perturb;

double tau_max, tau_min;
double *tmp_tau, *tmp_weight, *tmp_vel;

FILE *fcloud_out = NULL;
int icr_cloud_save = 1;

int postprc;

/* GSL */
const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;

gsl_interp_accel *gsl_acc;
gsl_interp  *gsl_linear;