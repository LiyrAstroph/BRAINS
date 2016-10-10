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

#include "dnestvars.h"
#include "allvars.h"
#include "proto.h"

#include "dnest_con.h"

void *best_model_thismodel, *best_model_std_thismodel;

void postprocess_con()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  double temperature = 1.0;
  double *pm, *pmstd;
  int num_ps, i, j;
  void *posterior_sample, *post_model;
  double *Fcon_mean;
  
  best_model_thismodel = malloc(size_of_modeltype);
  best_model_std_thismodel = malloc(size_of_modeltype);

  temperature = 1.0;
  dnest_postprocess(temperature);
  
  if(thistask == roottask)
  {
    FILE *fp, *fcon;
    char fname[200];
    // get number of lines in posterior sample file
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

        //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    //file for continuum reconstruction
    fcon = fopen("data/con_rec.txt", "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/con_rec.txt.\n");
      exit(0);
    }

    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);
    Fcon_mean = malloc(parset.n_con_recon * sizeof(double));

    for(j=0; j<parset.n_con_recon; j++)
      Fcon_mean[j] = 0.0;

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);

      calculate_con_from_model(post_model);

      if(gsl_rng_uniform(gsl_r) < 1.0)
      {
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%f %f\n", Tcon[j], Fcon[j]/con_scale);
        }
        fprintf(fcon, "\n");
      }

      for(j=0; j<parset.n_con_recon; j++)
      {
        Fcon_mean[j] += Fcon[j];
      }
    }

    fclose(fp);
    fclose(fcon);
    
    for(j=0; j<parset.n_con_recon; j++)
      Fcon_mean[j] /= num_ps;

    // mean continuum
    sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f\n", Tcon[i], Fcon_mean[i]/con_scale);
    }
    fclose(fp);  

    pm = (double *)best_model_thismodel;
    pmstd = (double *)best_model_std_thismodel;
    for(j=0; j<num_params; j++)
    {
      pm[j] = pmstd[j] = 0.0;
    }
    for(i=0; i<num_ps; i++)
    {
      for(j =0; j<num_params; j++)
        pm[j] += *((double *)posterior_sample + i*num_params + j );
    }

    for(j=0; j<num_params; j++)
      pm[j] /= num_ps;

    for(i=0; i<num_ps; i++)
    {
      for(j=0; j<num_params; j++)
        pmstd[j] += pow( *((double *)posterior_sample + i*num_params + j ) - pm[j], 2.0 );
    }

    for(j=0; j<num_params; j++)
    {
      if(num_ps > 1)
        pmstd[j] = sqrt(pmstd[j]/(num_ps-1.0));
      else
        pmstd[j] = 0.0;
    }  

    for(j = 0; j<num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_thismodel + j), *((double *)best_model_std_thismodel+j) ); 

    free(post_model);
    free(posterior_sample);
    free(Fcon_mean);
  }
  return;
}

void reconstruct_con()
{
  char *argv[]={""};

  reconstruct_con_init();
  dnest_con(0, argv);

  postprocess_con();

  /*if(thistask == roottask)
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
  }*/

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
  
  sigma = exp(pm[1]);
  tau = exp(pm[2]);
  alpha = 1.0;
  mu = pm[3];
  
  set_covar_Pmat(sigma, tau, alpha);
  Chol_decomp_L(PSmat, parset.n_con_recon, &info);
  multiply_matvec(PSmat, &pm[num_params_var], parset.n_con_recon, Fcon);

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
  double prob = 0.0, fcon, var2;
  int i;
  double *pm = (double *)model;
  
  calculate_con_from_model(model);
  
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  
  for(i=0; i<n_con_data; i++)
  {
    fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
    var2 = Fcerrs_data[i] * Fcerrs_data[i] + exp(pm[0])*exp(pm[0]);
    prob += -0.5*pow( (fcon - Fcon_data[i])/Fcerrs_data[i] ,  2.0) - 0.5*(log(2.0*PI) + log(var2));
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
  int i;
  double dT;

  /* set time array for continuum */
  Tcon_min = Tcon_data[0] - fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), 10.0);
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*(Tcon_data[n_con_data -1] - Tcon_data[0]), 10.0);
  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  Fcon = malloc(parset.n_con_recon * sizeof(double));

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSCON");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  perturb_accept = malloc(parset.num_particles * sizeof(int));
  for(i=0; i<parset.num_particles; i++)
    perturb_accept[i] = 0;

  return;
}

void reconstruct_con_end()
{
  free(Fcon);
  free(perturb_accept);
  return;
}