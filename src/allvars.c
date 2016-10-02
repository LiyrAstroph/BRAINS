/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
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


double VelUnit;

/* Data */
int n_con_data, n_line_data, n_vel_data;
double *Tcon_data, *Fcon_data,  *Fcerrs_data;
double *Tline_data, *Fline_data, *Flerrs_data;
double *Vline_data, *Fline2d_data, *Flerrs2d_data;
double con_scale, line_scale;

char dnest_options_file[BRAINS_MAX_STR_LENGTH];

int which_parameter_update, which_particle_update;  // which parameter and particle to be updated
int which_level_update;
int *perturb_accept, *which_parameter_update_prev;

/* continuum reconstruction */
double *Tcon, *Fcon, *Fcerrs;
double Tcon_min, Tcon_max;
double *PSmat, *USmat, *PSmat_data;

/* line reconstruction */
double *Fline_at_data;
double *Tline, *Fline, *Flerrs;
double Tline_min, Tline_max;

/* line reconstruction */
double *Fline2d_at_data;
double *Fline2d, *Flerrs2d;

// BLR
BLRmodel range_model[2];
int *par_fix, npar_fix;
double *par_fix_val;

// continuum variation
double var_range_model[4][2]; 

/* transfer function / velocity-delay map */
double *TransTau, *TransV, *Trans1D, *Trans2D_at_veldata, *Trans2D;
double dTransTau, dTransV;
double rcloud_min_set, rcloud_max_set;

double **Fcon_particles, **Fcon_particles_perturb;
double **Trans1D_particles, **Trans2D_at_veldata_particles;
double **Trans1D_particles_perturb, **Trans2D_at_veldata_particles_perturb;
double **clouds_particles, **clouds_particles_perturb;
double *prob_con_particles, *prob_con_particles_perturb;

/* GSL */
const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;

gsl_interp_accel *gsl_acc;
gsl_interp  *gsl_linear;