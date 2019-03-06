/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_con.c
 *  \brief reconstruct continuum.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>

#include "dnestvars.h"
#include "allvars.h"
#include "proto.h"

#include "dnest_con.h"

void *best_model_con;   /*!< best model */
void *best_model_std_con;  /*!< standard deviation of the best model */

/*!
 *  this function does postprocess for continuum. 
 */
void postprocess_con()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  double *pm, *pmstd;
  int num_ps, i, j;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);
  
  best_model_con = malloc(size_of_modeltype);
  best_model_std_con = malloc(size_of_modeltype);
  
  if(thistask == roottask)
  {
    char fname[200];
    FILE *fp, *fcon;

    /* get file name of posterior sample file */
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

    /* open file for posterior sample */
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }


    /* open file for continuum reconstruction */
    sprintf(fname, "%s/%s", parset.file_dir, "data/con_rec.txt");
    fcon = fopen(fname, "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }


    /* read number of points in posterior sample */
    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    which_parameter_update = -1;
    which_particle_update = 0;
    Fcon = Fcon_particles[which_particle_update];

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

      calculate_con_from_model_semiseparable(post_model);

      if(gsl_rng_uniform(gsl_r) < 1.0)
      {
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%f %f %f\n", Tcon[j], Fcon[j]/con_scale, Fcerrs[j]/con_scale);
        }
        fprintf(fcon, "\n");
      }
    }

    fclose(fp);
    fclose(fcon);

    /* calcaulte mean and standard deviation of posterior samples. */
    pm = (double *)best_model_con;
    pmstd = (double *)best_model_std_con;
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
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_con + j), 
                                             *((double *)best_model_std_con+j) ); 

    free(post_model);
    free(posterior_sample);
  }
  return;
}

/*! 
 * This function runs dnest sampling and postprocessing, reconstructs the continuum 
 *  using the best estimates of model parameters.
 */
void reconstruct_con()
{
  int i, argc=0;
  char **argv;

  // configure restart of dnest 
  argv = malloc(9*sizeof(char *));
  for(i=0; i<9; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restart_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart_dnest.txt");
  }
  if(parset.flag_postprc == 1)
  {
    strcpy(argv[argc++], "-p");
  }
  if(parset.flag_temp == 1)
  {
    sprintf(argv[argc++], "-t%f", parset.temperature);
  }
  if(parset.flag_sample_info == 1)
  {
    strcpy(argv[argc++], "-c");
  }
  
  //level-dependent sampling
  {
    strcpy(argv[argc++], "-l");
  }

  reconstruct_con_init();
  logz_con = dnest_con(argc, argv);

  if(parset.flag_exam_prior != 1)
  {
    postprocess_con();
  
    if(thistask == roottask)
    {
      which_parameter_update = -1;
      which_particle_update = 0;
      Fcon = Fcon_particles[which_particle_update];

      calculate_con_from_model_semiseparable(best_model_con);
 
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
        fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i] / con_scale, Fcerrs[i]/con_scale);
      }
      fclose(fp);
 
      memcpy(var_param, best_model_con, num_params_var*sizeof(double));
      memcpy(var_param_std, best_model_std_con, num_params_var*sizeof(double));
    }
  }
  else
  {
    for(i=0; i<num_params_var; i++)
    {
      var_param[i] = 0.5*(par_range_model[i][0] + par_range_model[i][1]);
      var_param_std[i] =0.5*(par_range_model[i][1] - par_range_model[i][0])/2.35;
    }
  }

  //use the posterior mean and standard variance as the prior for the 1d and 2d RM.
  //only for variability parameters.
  MPI_Bcast(var_param, num_params_var, MPI_DOUBLE, roottask, MPI_COMM_WORLD);
  MPI_Bcast(var_param_std, num_params_var, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  reconstruct_con_end();

  //clear up argv 
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

/*!
 *  this function calculates continuum light curves form model parameters 
 *  and stores it into Fcon.
 */
void calculate_con_from_model(const void *model)
{
  double *Larr, *ybuf, *y, *yu, *Cq, *yq, *yuq, *Larr_rec;
  double syserr;

  double *pm = (double *)model;
  double sigma, tau, alpha;
  int i, j, info;

  syserr = (exp(pm[0]) - 1.0) * con_error_mean;  // systematic error 
  tau = exp(pm[2]);
  sigma = exp(pm[1]) * sqrt(tau);
  alpha = 1.0;
  
  Larr = workspace; 
  ybuf = Larr + n_con_data*nq; 
  y = ybuf + n_con_data*nq;
  Cq = y + n_con_data;
  yq = Cq + nq*nq;
  yu = yq + nq; 
  yuq = yu + parset.n_con_recon;
  Larr_rec = yuq + parset.n_con_recon;

  set_covar_Pmat_data(sigma, tau, alpha, syserr);
  set_covar_Umat(sigma, tau, alpha);

  inverse_mat(PCmat_data, n_con_data, &info);

  for(i=0;i<n_con_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(Tcon_data[i], j);
  }
 
  // Cq^-1 = L^TxC^-1xL
  multiply_mat_MN(PCmat_data, Larr, ybuf, n_con_data, nq, n_con_data);
  multiply_mat_MN_transposeA(Larr, ybuf, Cq, nq, nq, n_con_data);

  // L^TxC^-1xy
  multiply_matvec(PCmat_data, Fcon_data, n_con_data, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

  // (hat q) = Cqx(L^TxC^-1xy)
  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  // q = uq + (hat q)
  Chol_decomp_L(Cq, nq, &info);
  multiply_matvec(Cq, &pm[3], nq, yq);
  for(i=0; i<nq; i++)
    yq[i] += ybuf[i];

  memcpy(con_q, yq, nq*sizeof(double)); //back up long-term trend
  
  // y = yc - Lxq
  multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
  for(i=0; i<n_con_data; i++)
  {
    y[i] = Fcon_data[i] - ybuf[i];
  }

  // (hat s) = SxC^-1xy
  multiply_matvec(PCmat_data, y, n_con_data, ybuf);
  multiply_matvec_MN(USmat, parset.n_con_recon, n_con_data, ybuf, Fcon);

  // SxC^-1xS^T
  multiply_mat_MN(USmat, PCmat_data, PEmat1, parset.n_con_recon, n_con_data, n_con_data);
  multiply_mat_MN_transposeB(PEmat1, USmat, PEmat2, parset.n_con_recon, parset.n_con_recon, n_con_data);

  set_covar_Pmat(sigma, tau, alpha);
  inverse_mat(PSmat, parset.n_con_recon, &info);

  // Q = [S^-1 + N^-1]^-1
  memcpy(PQmat, PSmat, parset.n_con_recon*parset.n_con_recon*sizeof(double));
  for(i=0; i<parset.n_con_recon; i++)
    PQmat[i*parset.n_con_recon+i] += 1.0/(sigma*sigma + syserr*syserr - PEmat2[i*parset.n_con_recon + i]);  
  
  inverse_mat(PQmat, parset.n_con_recon, &info);
  Chol_decomp_L(PQmat, parset.n_con_recon, &info);
  multiply_matvec(PQmat, &pm[num_params_var], parset.n_con_recon, yu);

  // add back long-term trend of continuum
  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }
  multiply_matvec_MN(Larr_rec, parset.n_con_recon, nq, yq, yuq);

  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcon[i] += yu[i] + yuq[i];
    Fcerrs[i] = sqrt(sigma*sigma + syserr*syserr - PEmat2[i*parset.n_con_recon + i]);
  }

  return;
}

/*!
 *  this function calculates continuum light curves form model parameters 
 *  and stores it into Fcon.
 *  
 *  using fast algorithm for semiseparable matrices
 */
void calculate_con_from_model_semiseparable(const void *model)
{
  double *Larr, *Lbuf, *ybuf, *y, *yu, *Cq, *yq, *yuq, *Larr_rec, *W, *D, *phi;
  double syserr;

  double *pm = (double *)model;
  double sigma, sigma2, tau, alpha;
  int i, j, info;

  syserr = (exp(pm[0]) - 1.0) * con_error_mean;  // systematic error 
  tau = exp(pm[2]);
  sigma = exp(pm[1]) * sqrt(tau);
  sigma2 = sigma*sigma;
  alpha = 1.0;
  
  Larr = workspace; 
  Lbuf = Larr + n_con_data*nq;
  ybuf = Lbuf + n_con_data*nq; 
  y = ybuf + n_con_data;
  Cq = y + n_con_data;
  yq = Cq + nq*nq;
  yu = yq + nq; 
  yuq = yu + parset.n_con_recon;
  Larr_rec = yuq + parset.n_con_recon;
  W = Larr_rec + parset.n_con_recon*nq;
  D = W + n_con_max;
  phi = D + n_con_max;

  for(i=0;i<n_con_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(Tcon_data[i], j);
  }
 
  set_covar_Umat(sigma, tau, alpha);

  compute_semiseparable_drw(Tcon_data, n_con_data, sigma2, 1.0/tau, Fcerrs_data, syserr, W, D, phi);
  // Cq^-1 = L^TxC^-1xL
  multiply_mat_semiseparable_drw(Larr, W, D, phi, n_con_data, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

  // L^TxC^-1xy
  multiply_matvec_semiseparable_drw(Fcon_data, W, D, phi, n_con_data, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

  // (hat q) = Cqx(L^TxC^-1xy)
  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  // q = uq + (hat q)
  Chol_decomp_L(Cq, nq, &info);
  multiply_matvec(Cq, &pm[3], nq, yq);
  for(i=0; i<nq; i++)
    yq[i] += ybuf[i];

  memcpy(con_q, yq, nq*sizeof(double)); //back up long-term trend
  
  // y = yc - Lxq
  multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
  for(i=0; i<n_con_data; i++)
  {
    y[i] = Fcon_data[i] - ybuf[i];
  }

  // (hat s) = SxC^-1xy
  multiply_matvec_semiseparable_drw(y, W, D, phi, n_con_data, sigma2, ybuf);
  multiply_matvec_MN(USmat, parset.n_con_recon, n_con_data, ybuf, Fcon);

  // SxC^-1xS^T
  multiply_mat_transposeB_semiseparable_drw(USmat, W, D, phi, n_con_data, parset.n_con_recon, sigma2, PEmat1);
  multiply_mat_MN(USmat, PEmat1, PEmat2, parset.n_con_recon, parset.n_con_recon, n_con_data);

  set_covar_Pmat(sigma, tau, alpha);

  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcerrs[i] = sqrt(sigma*sigma + syserr*syserr - PEmat2[i*parset.n_con_recon + i]);
  }
  compute_semiseparable_drw(Tcon, parset.n_con_recon, sigma2, 1.0/tau, Fcerrs, 0.0, W, D, phi);

  // Q = [S^-1 + N^-1]^-1 = N x [S+N]^-1 x S
  multiply_mat_semiseparable_drw(PSmat, W, D, phi, parset.n_con_recon, parset.n_con_recon, sigma2, PEmat2);
  for(i=0; i<parset.n_con_recon; i++)
  {
    for(j=0; j<=i; j++)
    {
      PQmat[i*parset.n_con_recon + j] = PQmat[j*parset.n_con_recon + i] = Fcerrs[i] * Fcerrs[i] * PEmat2[i*parset.n_con_recon+j];
    }
  }  

  Chol_decomp_L(PQmat, parset.n_con_recon, &info);
  multiply_matvec(PQmat, &pm[num_params_var], parset.n_con_recon, yu);

  // add back long-term trend of continuum
  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }
  multiply_matvec_MN(Larr_rec, parset.n_con_recon, nq, yq, yuq);

  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcon[i] += yu[i] + yuq[i];
  }

  return;
}


/*!
 * This function calculate continuum ligth curves from varibility parameters.
 */
void reconstruct_con_from_varmodel(double sigma_hat, double tau, double alpha, double syserr)
{
  double *Larr, *Lbuf, *ybuf, *y, *Larr_rec, *yq, *yuq, *Cq, sigma;
  int i, j, info;
  double *PEmat3, *PEmat4;

  sigma = sigma_hat * sqrt(tau);

  Larr = workspace;
  Lbuf = Larr + n_con_data * nq;
  ybuf = Lbuf + n_con_data * nq;
  y = ybuf + n_con_data;
  Cq = y + n_con_data;
  yq = Cq + nq*nq;
  yuq = yq + nq; 
  Larr_rec = yuq + parset.n_con_recon;

  PEmat3 = malloc(parset.n_con_recon * nq * sizeof(double));
  PEmat4 = malloc(parset.n_con_recon * parset.n_con_recon * sizeof(double));

  set_covar_Pmat_data(sigma, tau, alpha, syserr);
  set_covar_Umat(sigma, tau, alpha);

  inverse_mat(PCmat_data, n_con_data, &info);
  
  for(i=0;i<n_con_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(Tcon_data[i], j);
  }

  multiply_mat_MN(PCmat_data, Larr, Lbuf, n_con_data, nq, n_con_data);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

  multiply_matvec(PCmat_data, Fcon_data, n_con_data, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  memcpy(yq, ybuf, nq*sizeof(double));
  memcpy(con_q, ybuf, nq*sizeof(double)); // back up long-term trend
  
  multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
  for(i=0; i<n_con_data; i++)
  {
    y[i] = Fcon_data[i] - ybuf[i];
  }

  multiply_matvec(PCmat_data, y, n_con_data, ybuf);
  multiply_matvec_MN(USmat, parset.n_con_recon, n_con_data, ybuf, Fcon);

  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }
  multiply_matvec_MN(Larr_rec, parset.n_con_recon, nq, yq, yuq);

  for(i=0; i<parset.n_con_recon; i++)
    Fcon[i] += yuq[i];
  
  // SxC^-1xS^T
  multiply_mat_MN(USmat, PCmat_data, PEmat1, parset.n_con_recon, n_con_data, n_con_data);
  multiply_mat_MN_transposeB(PEmat1, USmat, PEmat2, parset.n_con_recon, parset.n_con_recon, n_con_data);

  multiply_mat_MN(PEmat1, Larr, PEmat3, parset.n_con_recon, nq, n_con_data);
  for(i=0; i<parset.n_con_recon*nq; i++)PEmat3[i] -= Larr_rec[i];
  multiply_mat_MN(PEmat3, Cq, PEmat1, parset.n_con_recon, nq, nq);
  multiply_mat_MN_transposeB(PEmat1, PEmat3, PEmat4, parset.n_con_recon, parset.n_con_recon, nq);

  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcerrs[i] = sqrt(sigma*sigma - PEmat2[i*parset.n_con_recon+i] + PEmat4[i*parset.n_con_recon + i]);
  }

  free(PEmat3);
  free(PEmat4);
  return;
}

/*!
 * this function calculates likelihood 
 */
double prob_con_variability(const void *model)
{
  double prob = 0.0;
  int i, j, param, info;
  double *pm = (double *)model;
  double tau, sigma, alpha, lndet, syserr;
  double *Larr, *Lbuf, *ybuf, *y, *yq, *Cq;
  
  which_particle_update = dnest_get_which_particle_update();
  
  if( which_parameter_update < num_params_var)
  {
    syserr = (exp(pm[0])-1.0)*con_error_mean;
    tau = exp(pm[2]);
    sigma = exp(pm[1]) * sqrt(tau);
    alpha = 1.0;
  
    Larr = workspace;
    Lbuf = Larr + n_con_data*nq;
    ybuf = Lbuf + n_con_data*nq;
    y = ybuf + n_con_data;
    yq = y + n_con_data;
    Cq = yq + nq;

    set_covar_Pmat_data(sigma, tau, alpha, syserr);
    memcpy(IPCmat_data, PCmat_data, n_con_data*n_con_data*sizeof(double));

    lndet = lndet_mat(PCmat_data, n_con_data, &info);

    inverse_mat(IPCmat_data, n_con_data, &info); /* calculate C^-1 */

    for(i=0;i<n_con_data;i++)
    {
      Larr[i*nq + 0]=1.0;
      for(j=1; j<nq; j++)
        Larr[i*nq + j] = pow(Tcon_data[i], j);
    }
 
    /* calculate L^T*C^-1*L */
    multiply_mat_MN(IPCmat_data, Larr, Lbuf, n_con_data, nq, n_con_data);
    multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

    /* calculate L^T*C^-1*y */
    multiply_matvec(IPCmat_data, Fcon_data, n_con_data, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_mat(Cq, nq, &info);
    multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

    Chol_decomp_L(Cq, nq, &info);
    multiply_matvec(Cq, &pm[3], nq, yq);
    for(i=0; i<nq; i++)
      yq[i] += ybuf[i];
  
    multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
    for(i=0; i<n_con_data; i++)
    {
      y[i] = Fcon_data[i] - ybuf[i];
    }
  
    multiply_matvec(IPCmat_data, y, n_con_data, ybuf);
    prob = -0.5 * cblas_ddot(n_con_data, y, 1, ybuf, 1);
    prob += -0.5*lndet;

    prob_con_particles_perturb[which_particle_update] = prob;
  }
  else
  {
    prob = prob_con_particles[which_particle_update];
  }

  return prob;
}

/*!
 * this function calculates likelihood 
 *  
 * using fast algorithm for semiseparable matrices
 */
double prob_con_variability_semiseparable(const void *model)
{
  double prob = 0.0;
  int i, j, param, info;
  double *pm = (double *)model;
  double tau, sigma2, alpha, lndet, syserr;
  double *Larr, *Lbuf, *ybuf, *y, *yq, *Cq, *W, *D, *phi;
  
  which_particle_update = dnest_get_which_particle_update();
  
  if( which_parameter_update < num_params_var)
  {

    syserr = (exp(pm[0])-1.0)*con_error_mean;
    tau = exp(pm[2]);
    sigma2 = exp(2.0*pm[1]) * tau;
    alpha = 1.0;
  
    Larr = workspace;
    Lbuf = Larr + n_con_data*nq;
    ybuf = Lbuf + n_con_data*nq;
    y = ybuf + n_con_data;
    yq = y + n_con_data;
    Cq = yq + nq;
    W = Cq + nq*nq;
    D = W + n_con_data;
    phi = D + n_con_data;

    for(i=0;i<n_con_data;i++)
    {
      Larr[i*nq + 0]=1.0;
      for(j=1; j<nq; j++)
        Larr[i*nq + j] = pow(Tcon_data[i], j);
    }
 
    compute_semiseparable_drw(Tcon_data, n_con_data, sigma2, 1.0/tau, Fcerrs_data, syserr, W, D, phi);
    lndet = 0.0;
    for(i=0; i<n_con_data; i++)
      lndet += log(D[i]);

    /* calculate L^T*C^-1*L */
    multiply_mat_semiseparable_drw(Larr, W, D, phi, n_con_data, nq, sigma2, Lbuf);
    multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

    /* calculate L^T*C^-1*y */
    multiply_matvec_semiseparable_drw(Fcon_data, W, D, phi, n_con_data, sigma2, ybuf);
    multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

    /* calculate (L^T*C^-1*L)^-1 * L^T*C^-1*y */
    inverse_mat(Cq, nq, &info);
    multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

    Chol_decomp_L(Cq, nq, &info);
    multiply_matvec(Cq, &pm[3], nq, yq);
    for(i=0; i<nq; i++)
      yq[i] += ybuf[i];
  
    multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
    for(i=0; i<n_con_data; i++)
    {
      y[i] = Fcon_data[i] - ybuf[i];
    }
  
    /* y^T x C^-1 x y*/
    multiply_matvec_semiseparable_drw(y, W, D, phi, n_con_data, sigma2, ybuf);
    prob = -0.5 * cblas_ddot(n_con_data, y, 1, ybuf, 1);
    prob += -0.5*lndet;

    prob_con_particles_perturb[which_particle_update] = prob;
  }
  else
  {
    prob = prob_con_particles[which_particle_update];
  }

  return prob;
}

/*!
 * this function calculates likelihood at initital step.
 */
double prob_con_variability_initial(const void *model)
{
  double prob = 0.0;
  int i, j, info;
  double *pm = (double *)model;
  double tau, sigma, alpha, lndet, syserr;
  double *Larr, *Lbuf, *ybuf, *y, *yq, *Cq;

  syserr = (exp(pm[0])-1.0)*con_error_mean;
  tau = exp(pm[2]);
  sigma = exp(pm[1]) * sqrt(tau);
  alpha = 1.0;
  
  Larr = workspace;
  Lbuf = Larr + n_con_data*nq;
  ybuf = Lbuf + n_con_data*nq;
  y = ybuf + n_con_data;
  yq = y + n_con_data;
  Cq = yq + nq;

  set_covar_Pmat_data(sigma, tau, alpha, syserr);
  memcpy(IPCmat_data, PCmat_data, n_con_data*n_con_data*sizeof(double));

  lndet = lndet_mat(PCmat_data, n_con_data, &info);

  inverse_mat(IPCmat_data, n_con_data, &info);

  for(i=0;i<n_con_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(Tcon_data[i], j);
  }
 
  multiply_mat_MN(IPCmat_data, Larr, Lbuf, n_con_data, nq, n_con_data);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

  multiply_matvec(IPCmat_data, Fcon_data, n_con_data, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  Chol_decomp_L(Cq, nq, &info);
  multiply_matvec(Cq, &pm[3], nq, yq);
  for(i=0; i<nq; i++)
    yq[i] += ybuf[i];
  
  multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
  for(i=0; i<n_con_data; i++)
  {
    y[i] = Fcon_data[i] - ybuf[i];
  }
  
  multiply_matvec(IPCmat_data, y, n_con_data, ybuf);
  prob = -0.5 * cblas_ddot(n_con_data, y, 1, ybuf, 1);
  prob += -0.5*lndet;

  prob_con_particles[which_particle_update] = prob;
  return prob;
}

/*!
 * this function calculates likelihood at initital step.
 *  
 * using fast algorithm for semiseparable matrices
 */
double prob_con_variability_initial_semiseparable(const void *model)
{
  double prob = 0.0;
  int i, j, info;
  double *pm = (double *)model;
  double tau, sigma2, alpha, lndet, syserr;
  double *Larr, *Lbuf, *ybuf, *y, *yq, *Cq, *W, *D, *phi;

  syserr = (exp(pm[0])-1.0)*con_error_mean;
  tau = exp(pm[2]);
  sigma2 = exp(2.0*pm[1]) * tau;
  alpha = 1.0;
  
  Larr = workspace;
  Lbuf = Larr + n_con_data*nq;
  ybuf = Lbuf + n_con_data*nq;
  y = ybuf + n_con_data;
  yq = y + n_con_data;
  Cq = yq + nq;
  W = Cq + nq*nq;
  D = W + n_con_data;
  phi = D + n_con_data;

  for(i=0;i<n_con_data;i++)
  {
    Larr[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr[i*nq + j] = pow(Tcon_data[i], j);
  }

  compute_semiseparable_drw(Tcon_data, n_con_data, sigma2, 1.0/tau, Fcerrs_data, syserr, W, D, phi);
  lndet = 0.0;
  for(i=0; i<n_con_data; i++)
  lndet += log(D[i]);
 
  /* calculate L^T*C^-1*L */
  multiply_mat_semiseparable_drw(Larr, W, D, phi, n_con_data, nq, sigma2, Lbuf);
  multiply_mat_MN_transposeA(Larr, Lbuf, Cq, nq, nq, n_con_data);

  /* calculate L^T*C^-1*y */
  multiply_matvec_semiseparable_drw(Fcon_data, W, D, phi, n_con_data, sigma2, ybuf);
  multiply_mat_MN_transposeA(Larr, ybuf, yq, nq, 1, n_con_data);

  inverse_mat(Cq, nq, &info);
  multiply_mat_MN(Cq, yq, ybuf, nq, 1, nq);

  Chol_decomp_L(Cq, nq, &info);
  multiply_matvec(Cq, &pm[3], nq, yq);
  for(i=0; i<nq; i++)
    yq[i] += ybuf[i];
  
  multiply_matvec_MN(Larr, n_con_data, nq, yq, ybuf);
  for(i=0; i<n_con_data; i++)
  {
    y[i] = Fcon_data[i] - ybuf[i];
  }
  
  /* y^T x C^-1 x y*/
  multiply_matvec_semiseparable_drw(y, W, D, phi, n_con_data, sigma2, ybuf);
  prob = -0.5 * cblas_ddot(n_con_data, y, 1, ybuf, 1);
  prob += -0.5*lndet;

  prob_con_particles[which_particle_update] = prob;
  return prob;
}

/*!
 * this function sets the covariance matrix at time points for reconstruction.
 * store into PSmat.
 */
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

/*!
 * this function sets the covariance matrix at data time points 
 */
void set_covar_Pmat_data(double sigma, double tau, double alpha, double syserr)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<n_con_data; i++)
  { 
    t1 = Tcon_data[i];

    for(j=0; j<i; j++)
    {
      t2 = Tcon_data[j];
      PSmat_data[i*n_con_data+j] = sigma*sigma* exp (- pow (fabs(t1-t2) / tau, alpha));
      PSmat_data[j*n_con_data+i] = PSmat_data[i*n_con_data+j];

      PNmat_data[i*n_con_data+j] = PNmat_data[j*n_con_data+i] = 0.0;

      PCmat_data[i*n_con_data+j] = PCmat_data[j*n_con_data+i] = PSmat_data[i*n_con_data+j];
    }

    PSmat_data[i*n_con_data+i] = sigma * sigma;
    PNmat_data[i*n_con_data+i] = Fcerrs_data[i]*Fcerrs_data[i] + syserr*syserr;
    PCmat_data[i*n_con_data+i] = PSmat_data[i*n_con_data+i] + PNmat_data[i*n_con_data+i];
  }
  return;
}

/*!
 * this function sets the covariance matrix at time of data points and reconstruction points
 */
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

/*!
 * this function initializes the continuum reconstruction.
 */
void reconstruct_con_init()
{
  int i;
  double dT, Tspan;

  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];

  /* set time array for continuum */
  if(parset.time_back > 0.0)
    Tcon_min = Tcon_data[0] - parset.time_back; 
  else
    Tcon_min = Tcon_data[0] - fmax(0.05*Tspan, 10.0);

  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 10.0);

  if(thistask == 0)
    printf("Tcon_min_max: %f %f\n", Tcon_min - Tcon_data[0], Tcon_max - Tcon_data[n_con_data-1]);

  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSCON");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  Fcon_particles = malloc(parset.num_particles * sizeof(double *));
  Fcon_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fcon_particles[i] = malloc(parset.n_con_recon * sizeof(double));
    Fcon_particles_perturb[i] = malloc(parset.n_con_recon * sizeof(double));
  }
  prob_con_particles = malloc(parset.num_particles * sizeof(double));
  prob_con_particles_perturb = malloc(parset.num_particles * sizeof(double));

  con_q = malloc(nq * sizeof(double));
  return;
}

/*!
 * this function finalize the continuum reconstruction.
 */
void reconstruct_con_end()
{
  int i;

  for(i=0; i<parset.num_particles; i++)
  {
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
  }
  free(Fcon_particles);
  free(Fcon_particles_perturb);
  
  free(prob_con_particles);
  free(prob_con_particles_perturb);

  free(best_model_con);
  free(best_model_std_con);

  free(con_q);
  return;
}

/*!
 *  create continumm from random number generate based on damped random walk model.
 */

void create_con_from_random(double sigma_hat, double tau, double alpha, double syserr)
{
  int i, info;
  double *Prandvec;
  double sigma = sigma_hat * sqrt(tau);
  double mean;

  Prandvec = malloc(parset.n_con_recon*sizeof(double));

  set_covar_Pmat(sigma, tau, alpha);

  for(i=0; i<parset.n_con_recon; i++)
  {
    Prandvec[i] = gsl_ran_ugaussian(gsl_r);
  }

  Chol_decomp_L(PSmat, parset.n_con_recon, &info);

  multiply_matvec(PSmat, Prandvec, parset.n_con_recon, Fcon);
  
  mean = 0.0;
  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcon[i] += sigma*3.0; /* shift light curve to have positive fluxes*/
    mean += Fcon[i];
  }
  mean /= parset.n_con_recon;

  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcon[i] /= mean;
  }
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Fcon[i] +=  gsl_ran_ugaussian(gsl_r) * 0.01;
    Fcerrs[i] = 0.01;
  }

  free(Prandvec);
  return;
}