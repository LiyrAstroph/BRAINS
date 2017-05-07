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


double VelUnit;    /*!< Velocity unit = \f$\sqrt{\frac{G\cdot 10^6M_\odot}{1~\rm light-day}}/10^5\f$km/s, 
                        with \f$ M_\bullet = \f$ and \f$r=1\f$light day.*/

/* Data */
int n_con_data;    /*!< number of continuum data points */
int n_line_data;   /*!< number of emission data points along time axis */
int n_vel_data;    /*!< number of emission data points along velocity axis */
double *Tcon_data, *Fcon_data,  *Fcerrs_data;
double *Tline_data, *Fline_data, *Flerrs_data;
double *Vline_data, *Fline2d_data, *Flerrs2d_data;
double con_scale, line_scale;

char dnest_options_file[BRAINS_MAX_STR_LENGTH];

int which_parameter_update;  /*!< which parameter being updated */
int which_particle_update;   /*!< which particle being updated */
int which_level_update;      /*!< which level of the particle */
unsigned long long int which_mcmc_steps;//mcmc steps 
int *perturb_accept, *which_parameter_update_prev;

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


/* line reconstruction */
double *Fline_at_data;
double *Tline, *Fline, *Flerrs;
double Tline_min, Tline_max;

/* line reconstruction */
double *Fline2d_at_data;
double *Fline2d, *Flerrs2d;

// BLR
int BLRmodel_size;
int *par_fix, npar_fix;
double *par_fix_val;
double **blr_range_model, **par_range_model;

int num_params, num_params_blr, num_params_var, num_params_trend;

// continuum variation
double var_range_model[7][2]; /*!< define the range of variability parameters */

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
double *prob_line_particles, *prob_line_particles_perturb;

int postprc;

/* GSL */
const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;

gsl_interp_accel *gsl_acc;
gsl_interp  *gsl_linear;