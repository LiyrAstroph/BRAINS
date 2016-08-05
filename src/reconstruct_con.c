/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include "dnest_con.h"

void reconstruct_con()
{
  double *pm;
  double sigma, tau, alpha, mu;

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

}

/* calculate continuum light curves form model parameters 
 * stored into Fcon.
 */
void calculate_con_from_model(const void *model)
{
  double *pm = (double *)model;
  double prob = 0.0, sigma, tau, alpha, mu, fcon;
  int i, info;
  
  sigma = exp(pm[0]);
  tau = exp(pm[1]);
  alpha = 1.0;
  mu = pm[2];
  
  set_covar_Pmat(sigma, tau, alpha);
  Chol_decomp_U(PSmat, parset.n_con_recon, &info);
  multiply_matvec(PSmat, &pm[3], parset.n_con_recon, Fcon);

  // add back the mean of continuum
  for(i=0; i<parset.n_con_recon; i++)
    Fcon[i] += mu;

  return;
}

/* calculate likelihood */
double prob_con_variability(const void *model)
{
  double prob = 0.0, fcon;
  int i, info;
  
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

void reconstruct_con_init()
{
  return;
}