/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_interp.h>

#include "allvars.h"
#include "proto.h"

#include "dnest_con.h"

void reconstruct_con()
{
  reconstruct_con_init();
  dnest_con(0, " ");

  if(thistask == roottask)
  {
    calculate_con_from_model(best_model_thismodel);

    FILE *fp;
    char fname[200];
    int i;
    sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }

    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tcon[i], Fcon[i] / con_scale);
    }
    fclose(fp);
  }

  reconstruct_con_end();

}

/* calculate continuum light curves form model parameters 
 * stored into Fcon.
 */
void calculate_con_from_model(const void *model)
{
  double *pm = (double *)model;
  double sigma, tau, alpha, mu;
  int i, info;
  
  sigma = exp(pm[0]);
  tau = exp(pm[1]);
  alpha = 1.0;
  mu = pm[2];
  
  set_covar_Pmat(sigma, tau, alpha);
  Chol_decomp_L(PSmat, parset.n_con_recon, &info);
  multiply_matvec(PSmat, &pm[3], parset.n_con_recon, Fcon);

  // add back the mean of continuum
  for(i=0; i<parset.n_con_recon; i++)
    Fcon[i] += mu;

  return;
}

void reconstruct_con_from_varmodel(double sigma, double tau, double alpha)
{
  double *Larr, *ybuf, *y;
  double lambda, ave_con;
  int i, info;

  Larr = malloc(n_con_data * sizeof(double));
  ybuf = malloc(n_con_data* sizeof(double));
  y = malloc(n_con_data* sizeof(double));

  set_covar_Pmat_data(sigma, tau, alpha);
  set_covar_Umat(sigma, tau, alpha);

  inverse_mat(PSmat_data, n_con_data, &info);

  for(i=0;i<n_con_data;i++)
    Larr[i]=1.0;
  multiply_matvec(PSmat_data, Larr, n_con_data, ybuf);
  lambda = cblas_ddot(n_con_data, Larr, 1, ybuf, 1);
  multiply_matvec(PSmat_data, Fcon_data, n_con_data, ybuf);
  ave_con = cblas_ddot(n_con_data, Larr, 1, ybuf, 1);
  ave_con /= lambda;
  
  for(i=0; i<n_con_data; i++)
    y[i] = Fcon_data[i] - ave_con;
  multiply_matvec(PSmat_data, y, n_con_data, ybuf);
  multiply_matvec_MN(USmat, parset.n_con_recon, n_con_data, ybuf, Fcon);

  for(i=0; i<parset.n_con_recon; i++)
    Fcon[i] += ave_con;

  free(Larr);
  free(y);
  free(ybuf);
}

/* calculate likelihood */
double prob_con_variability(const void *model)
{
  double prob = 0.0, fcon;
  int i;
  
  calculate_con_from_model(model);
  
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  
  for(i=0; i<n_con_data; i++)
  {
    fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
    prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - ( 0.5*log(2.0*PI) + log(Fcerrs_data[i]) );
  }

  return prob;
}

/* set the covariance matrix */
void set_covar_Pmat(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<parset.n_con_recon; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<=i; j++)
    {
      t2 = Tcon[j];
      PSmat[i*parset.n_con_recon+j] = sigma*sigma* exp (- pow (fabs(t1-t2) / tau, alpha));
      PSmat[j*parset.n_con_recon+i] = PSmat[i*parset.n_con_recon+j];
    }
  }
  return;
}

/* set the covariance matrix at data points */
void set_covar_Pmat_data(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n_con_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<=i; j++)
    {
      t2 = Tcon_data[j];
      PSmat_data[i*n_con_data+j] = sigma*sigma* exp (- pow (fabs(t1-t2) / tau, alpha));
      PSmat_data[j*n_con_data+i] = PSmat_data[i*n_con_data+j];
    }
  }
  return;
}

void set_covar_Umat(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<parset.n_con_recon; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<n_con_data; j++)
    {
      t2 = Tcon_data[j];
      USmat[i*n_con_data+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}

void reconstruct_con_init()
{
  Fcon = malloc(parset.n_con_recon * sizeof(double));
  return;
}

void reconstruct_con_end()
{
  free(Fcon);
  return;
}