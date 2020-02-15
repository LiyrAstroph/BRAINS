/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_line2d.c
 *  \brief reconstruct 2d line and BLR model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <gsl/gsl_interp.h>

#include "brains.h"

void *best_model_line2d;      /*!< best model */
void *best_model_std_line2d;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess2d()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  int num_ps, i, j, k, nc;
  double *pm, *pmstd;
  double *lag;
  void *posterior_sample, *post_model;
  double mean_lag, mean_lag_std, sum1, sum2;
  int size_of_modeltype = num_params * sizeof(double);
  
  best_model_line2d = malloc(size_of_modeltype);
  best_model_std_line2d = malloc(size_of_modeltype);

  if(thistask == roottask)
  {
    // initialize smoothing workspace
    smooth_init(n_vel_data_ext, Vline_data_ext);
    char fname[200];
    FILE *fp, *fcon, *fline, *ftran, *fline1d;
    double *Fline1d, dV;

    // velocity grid width, in term of wavelength of Hbeta.
    dV = (Vline_data[n_vel_data-1]-Vline_data[0])/(n_vel_data-1) * VelUnit * 4861.0/3.0e5; 
    Fline1d = malloc(n_line_data * sizeof(double));

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
    sprintf(fname, "%s/%s", parset.file_dir, "data/con_rec.txt");
    fcon = fopen(fname, "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    //file for line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/line2d_rec.txt");
    fline = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/line_rec.txt");
    fline1d = fopen(fname, "w");
    if(fline1d == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran2d_rec.txt");
    ftran = fopen(fname, "w");
    if(ftran == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    // read number of lines
    if(fscanf(fp, "# %d", &num_ps) < 1)
    {
      fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
      exit(0);
    }
    printf("# Number of points in posterior sample: %d\n", num_ps);

    lag = malloc(num_ps * sizeof(double));
    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    force_update = 1;
    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;
    mean_lag = 0.0;
    nc = 0;

    con_q = con_q_particles[which_particle_update];
    Fcon = Fcon_particles[which_particle_update];
    TransTau = TransTau_particles[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
    Fline2d_at_data = Fline_at_data_particles[which_particle_update];

    for(i=0; i<num_ps; i++)
    {
      // read lines
      for(j=0; j<num_params; j++)
      {
        if(fscanf(fp, "%lf", (double *)post_model + j) < 1)
        {
          fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file);
          exit(0);
        }
      }
      fscanf(fp, "\n");

      //store model
      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);
    
      //calculate_con_from_model(post_model + num_params_blr *sizeof(double));
      calculate_con_from_model_semiseparable(post_model + num_params_blr *sizeof(double));
      gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

      transfun_2d_cal(post_model, Vline_data_ext, Trans2D_at_veldata, 
                                      n_vel_data_ext, 0);
      calculate_line2d_from_blrmodel(post_model, Tline_data, Vline_data_ext, Trans2D_at_veldata, 
                                      Fline2d_at_data, n_line_data, n_vel_data_ext);
      
      // calculate integrated line fluxes
      for(j = 0; j < n_line_data; j++)
      {
        Fline1d[j] = 0.0;
        for(k=0; k<n_vel_data; k++)
        {
          Fline1d[j] += Fline2d_at_data[j * n_vel_data_ext + (k+n_vel_data_incr)] * dV /line_scale;
        }
      }

      // calculate mean time lag
      sum1 = 0.0;
      sum2 = 0.0;
      for(j=0; j<parset.n_tau; j++)
      {
        for(k=0; k<n_vel_data_ext; k++)
        {
          sum1 += Trans2D_at_veldata[j * n_vel_data_ext + k] * TransTau[j];
          sum2 += Trans2D_at_veldata[j * n_vel_data_ext + k];
        }
      }
      // take care of zero transfer function
      if(sum2 > 0.0)
      {
        lag[i] = sum1/sum2;
        mean_lag += lag[i];
        nc++;
      }
      else
      {
        lag[i] = -DBL_MAX;
      } 

      //if( i % (num_ps/10+1) == 0)  
      {
        // output continuum
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%e %e\n", Tcon[j]*(1.0+parset.redshift), Fcon[j]/con_scale);
        }
        fprintf(fcon, "\n");

        // output 2d line
        for(j=0; j<n_line_data; j++)
        {
          for(k=0; k<n_vel_data; k++)
          {
            fprintf(fline, "%e ", Fline2d_at_data[j * n_vel_data_ext + (k+n_vel_data_incr)]/line_scale);
          }
          fprintf(fline, "\n");
        }
        fprintf(fline, "\n");

        // output transfer function
        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran, "%e ", TransTau[j]);
          
          for(k=0; k<n_vel_data; k++)
          {
            fprintf(ftran, "%e ", Trans2D_at_veldata[j * n_vel_data_ext + (k+n_vel_data_incr)]);
          }
          fprintf(ftran, "\n");
        }
        fprintf(ftran, "\n");

        // output 1d line 
        for(j = 0; j<n_line_data; j++)
        {
          fprintf(fline1d, "%e %e\n", Tline_data[j]*(1.0+parset.redshift), Fline1d[j]);
        }
        fprintf(fline1d, "\n");
      }
    }

    smooth_end();

    fclose(fp);
    fclose(fcon);
    fclose(fline);
    fclose(ftran);
    fclose(fline1d);

    /* calculate mean lags */
    mean_lag /= (nc);
    mean_lag_std = 0.0;
    for(i=0; i<num_ps; i++)
    {
      if(lag[i] > -DBL_MAX)
        mean_lag_std += (lag[i] - mean_lag) * (lag[i] - mean_lag);
    }
    if(nc > 1)
      mean_lag_std = sqrt(mean_lag_std/(nc -1.0));
    else
      mean_lag_std = 0.0;
    printf("Mean time lag: %f+-%f\n", mean_lag, mean_lag_std);

    pm = (double *)best_model_line2d;
    pmstd = (double *)best_model_std_line2d;
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

    for(j = 0; j<num_params_blr + num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_line2d + j), 
                                             *((double *)best_model_std_line2d+j) ); 
 
    free(lag);
    free(post_model);
    free(posterior_sample);
    free(Fline1d);
  }
  return;
}

/*!
 * this function run dnest sampleing, reconstruct light curves using the best estimates for parameters.
 */
void reconstruct_line2d()
{
  int i, argc=0;
  char **argv;

  //configure restart of dnest
  argv = malloc(9*sizeof(char *));
  for(i=0; i<9; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restart2d_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restart2d_dnest.txt");
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
  
  reconstruct_line2d_init();
  
  smooth_init(n_vel_data_ext, Vline_data_ext);
  dnest_line2d(argc, argv);
  smooth_end();

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess2d();

    // calculate light curves using the best model
    if(thistask == roottask)
    {
      force_update = 1; 
      which_parameter_update = -1;
      which_particle_update = 0;
      con_q = con_q_particles[which_particle_update];
      Fcon = Fcon_particles[which_particle_update];
      TransTau = TransTau_particles[which_particle_update];
      Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
      Fline2d_at_data = Fline_at_data_particles[which_particle_update];
    
      //calculate_con_from_model(best_model_line2d + num_params_blr *sizeof(double));
      calculate_con_from_model_semiseparable(best_model_line2d + num_params_blr *sizeof(double));
      gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

      FILE *fp;
      char fname[200];
      int i, j;

      sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      for(i=0; i<parset.n_con_recon; i++)
      {
        fprintf(fp, "%e %e\n", Tcon[i]*(1.0+parset.redshift), Fcon[i] / con_scale);
      }
      fclose(fp);

      smooth_init(n_vel_data_ext, Vline_data_ext);
      // recovered line2d at data points
      transfun_2d_cal(best_model_line2d, Vline_data_ext, Trans2D_at_veldata, 
                                              n_vel_data_ext, parset.flag_save_clouds);
      calculate_line2d_from_blrmodel(best_model_line2d, Tline_data, Vline_data_ext, Trans2D_at_veldata, 
                                                       Fline2d_at_data, n_line_data, n_vel_data_ext);
    
      sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_data_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
    
      for(i=0; i<n_line_data; i++)
      {
        for(j=0; j<n_vel_data; j++)
        {
          fprintf(fp, "%e %e %e\n", Vline_data[j]*VelUnit, Tline_data[i]*(1.0+parset.redshift),  
                            Fline2d_at_data[i*n_vel_data_ext + (j+n_vel_data_incr)] / line_scale);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      sprintf(fname, "%s/%s", parset.file_dir, parset.tran2d_data_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
    
      fprintf(fp, "# %d %d\n", parset.n_tau, n_vel_data);
      for(i=0; i<parset.n_tau; i++)
      {
        for(j=0; j<n_vel_data; j++)
        {
          fprintf(fp, "%e %e %e\n", Vline_data[j]*VelUnit, TransTau[i],  Trans2D_at_veldata[i*n_vel_data_ext + (j+n_vel_data_incr)]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      smooth_end();
    
      // recovered line2d at specified points
      smooth_init(parset.n_vel_recon, TransV);

      which_parameter_update = -1;
      which_particle_update = 0;
      transfun_2d_cal(best_model_line2d, TransV, Trans2D, parset.n_vel_recon, 0);

      /* there is no data for spectral broadening at given specified epoch, using the mean value
       * and set InstRes_err=0.0.
       */ 
      double *pm = (double *)best_model_line2d;
      if(parset.flag_InstRes > 0) 
      {
        parset.flag_InstRes = 0; /* force to be uniform prior */
        double instres_mean = 0.0;
        for(i=0; i<n_line_data; i++)
        {
          instres_mean += instres_epoch[i];
        }
        parset.InstRes = instres_mean/n_line_data;
        parset.InstRes_err = 0.0;
        
        instres_mean = 0.0;
        for(i=0; i<n_line_data; i++)
        {
          instres_mean += pm[num_params_blr_model+num_params_nlr+i];
        }
        pm[num_params_blr_model + num_params_nlr ] = instres_mean/n_line_data;
      }

      /* similarly, there is no data for line center information at given specified epoch,
       * using the mean value
       */
      if(parset.flag_linecenter < 0)
      {
        parset.flag_linecenter = 1; /* force to be uniform prior, note num_params_linecenter is still unchanged */
        double linecenter_mean = 0.0;
        for(i=0; i<n_line_data; i++)
        {
          linecenter_mean += pm[num_params_blr - num_params_linecenter - 3 + i];
        }
        pm[num_params_blr - num_params_linecenter - 3] = linecenter_mean/n_line_data;
      }
      
      calculate_line2d_from_blrmodel(best_model_line2d, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

      sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      for(i=0; i<parset.n_line_recon; i++)
      {
        for(j=0; j<parset.n_vel_recon; j++)
        {
          fprintf(fp, "%e %e %e\n", TransV[j]*VelUnit, Tline[i],  Fline2d[i*parset.n_vel_recon + j] / line_scale);
        }

        fprintf(fp, "\n");
      }
      fclose(fp);

      // output 2d transfer function
      sprintf(fname, "%s/%s", parset.file_dir, parset.tran2d_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      fprintf(fp, "# %d %d\n", parset.n_tau, parset.n_vel_recon);
      for(i=0; i<parset.n_tau; i++)
      {
        for(j=0; j<parset.n_vel_recon; j++)
        {
          fprintf(fp, "%e %e %e\n", TransV[j]*VelUnit, TransTau[i], Trans2D[i*parset.n_vel_recon + j]);
        }

        fprintf(fp, "\n");
      }
      fclose(fp);
      smooth_end();
    }
  }

  reconstruct_line2d_end();

  //clear up argv
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

/*!
 * this function initializes 2d reconstruction.
 */
void reconstruct_line2d_init()
{
  int i, j;
  double dT, Tspan;

  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];

  /* set time grid for continuum */
  /*if(parset.time_back > 0.0)
    Tcon_min = Tcon_data[0] - parset.time_back;
  else
    Tcon_min = Tcon_data[0] - fmax(0.05*Tspan, Tspan/2.0 + (Tcon_data[0] - Tline_data[0]));*/
  
  Tcon_min = Tcon_data[0] - time_back_set;

  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 10.0);
  Tcon_max = fmax(Tcon_max, Tline_data[n_line_data -1]);  /* The time span should cover that of the emission line data */
  
  if(thistask == roottask)
    printf("Tcon_min_max: %f %f\n", Tcon_min - Tcon_data[0], Tcon_max - Tcon_data[n_con_data-1]);

  dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    Tcon[i] = Tcon_min + i*dT;
  }

  /* set Larr_rec */
  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }

  //TransTau = malloc(parset.n_tau * sizeof(double));
  //Trans2D_at_veldata = malloc(parset.n_tau * n_vel_data * sizeof(double));


  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  
  //Fline2d_at_data = malloc(n_line_data * n_vel_data * sizeof(double));

  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  
  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  if(parset.time_back <= 0.0)
    Tline_min = fmax(Tline_min, Tcon_min + Tspan/2.0 - 10.0);

  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10);
  Tline_max = fmin(Tline_max, Tcon_max);  /* The time span should be smaller than that of the continuum */

  if(thistask == roottask)
    printf("Tline_min_max: %f %f\n", Tline_min - Tline_data[0], Tline_max - Tline_data[n_line_data -1]);

  dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }

  double vel_max_set = Vline_data_ext[n_vel_data_ext -1], vel_min_set = Vline_data_ext[0];
  double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

  for(i=0; i<parset.n_vel_recon; i++)
  {
  	TransV[i] = vel_min_set + dVel*i;
  }

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONS2D");
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

  TransTau_particles = malloc(parset.num_particles * sizeof(double *));
  TransTau_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    TransTau_particles[i] = malloc(parset.n_tau * sizeof(double));
    TransTau_particles_perturb[i] = malloc(parset.n_tau * sizeof(double));
  }

  Trans2D_at_veldata_particles = malloc(parset.num_particles * sizeof(double *));
  Trans2D_at_veldata_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans2D_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_data_ext * sizeof(double));
    Trans2D_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_data_ext * sizeof(double));
  }

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_line_data * n_vel_data_ext * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_line_data * n_vel_data_ext * sizeof(double));
  }

  con_q_particles = malloc(parset.num_particles * sizeof(double *));
  con_q_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    con_q_particles[i] = malloc(nq * sizeof(double));
    con_q_particles_perturb[i] = malloc(nq* sizeof(double));
  }

  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));

  if(parset.flag_save_clouds && thistask == roottask)
  {
    if(parset.n_cloud_per_task <= 10000)
      icr_cloud_save = 1;
    else
      icr_cloud_save = parset.n_cloud_per_task/10000;

    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }
  return;
}

/*!
 * this function finalizes 2d reconstruction.
 */
void reconstruct_line2d_end()
{
  free(Tline);
  //free(Fline2d_at_data);
  free(Fline2d);

  //free(TransTau);
  free(TransV);
  free(Trans2D);
  //free(Trans2D_at_veldata);
  int i;
  for(i=0; i<parset.num_particles; i++)
  {
    free(Fcon_particles[i]);
    free(Fcon_particles_perturb[i]);
  }
  free(Fcon_particles);
  free(Fcon_particles_perturb);
  
  free(Fline_at_data);

  free(par_fix);
  free(par_fix_val);

  free(best_model_line2d);
  free(best_model_std_line2d);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans2D_at_veldata_particles[i]);
    free(Trans2D_at_veldata_particles_perturb[i]);
    free(TransTau_particles[i]);
    free(TransTau_particles_perturb[i]);
    free(Fline_at_data_particles[i]);
    free(Fline_at_data_particles_perturb[i]);

    free(con_q_particles[i]);
    free(con_q_particles_perturb[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(Trans2D_at_veldata_particles_perturb);
  free(TransTau_particles);
  free(TransTau_particles_perturb);
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);

  free(con_q_particles);
  free(con_q_particles_perturb);
 
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
    free(par_prior_gaussian[i]);
  }
  free(par_range_model);
  free(par_prior_gaussian);
  free(par_prior_model);
  
  free(clouds_tau);
  free(clouds_weight);
  free(clouds_vel);
  
  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  if(thistask == roottask)
  {
    printf("Ends reconstruct_line2d.\n");
  }
  return;
}

/*!
 * this function calculate probability at initial step.
 */
double prob_initial_line2d(const void *model)
{
  double prob_line = 0.0, var2, dy, var2_se;
  int i, j;
  double *pm = (double *)model;
  
  which_particle_update = dnest_get_which_particle_update();

  con_q = con_q_particles[which_particle_update];
  Fcon = Fcon_particles[which_particle_update];
  //calculate_con_from_model(model + num_params_blr*sizeof(double));
  calculate_con_from_model_semiseparable(model + num_params_blr*sizeof(double));
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

  TransTau = TransTau_particles[which_particle_update];
  Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
  Fline2d_at_data = Fline_at_data_particles[which_particle_update];
  which_parameter_update = -1;
  transfun_2d_cal(model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
  calculate_line2d_from_blrmodel(model, Tline_data, Vline_data_ext, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data_ext);

  var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean*line_error_mean;
  for(i=0; i<n_line_data; i++)
  {
    for(j=0; j<n_vel_data; j++)
    {
      //note mask with error < 0.0
      if(Flerrs2d_data[i] > 0.0)
      {
        dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
        var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }

  return prob_line;
}

/*!
 * this function calculate probability at restart step.
 */
double prob_restart_line2d(const void *model)
{
  double prob_line = 0.0, var2, dy, var2_se;
  int i, j;
  double *pm = (double *)model;
  
  which_particle_update = dnest_get_which_particle_update();
  con_q = con_q_particles[which_particle_update];
  Fcon = Fcon_particles[which_particle_update];
  //calculate_con_from_model(model + num_params_blr*sizeof(double));
  calculate_con_from_model_semiseparable(model + num_params_blr*sizeof(double));
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);

  TransTau = TransTau_particles[which_particle_update];
  Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
  Fline2d_at_data = Fline_at_data_particles[which_particle_update];
  which_parameter_update = num_params + 1; // so as not to update clouds.
  transfun_2d_cal(model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
  calculate_line2d_from_blrmodel(model, Tline_data, Vline_data_ext, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data_ext);

  var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean*line_error_mean;
  for(i=0; i<n_line_data; i++)
  {
    for(j=0; j<n_vel_data; j++)
    {
      //note mask with error < 0.0
      if(Flerrs2d_data[i] > 0.0)
      {
        dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
        var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }

  return prob_line;
}
/*!
 * this function calculate probability.
 * 
 * At each MCMC step, only one parameter is updated, which only changes some values; thus,
 * optimization that reuses the unchanged values can improve computation efficiency.
 */
double prob_line2d(const void *model)
{
  double prob_line = 0.0, var2, dy, var2_se;
  int i, j;
  double *pm = (double *)model;
  
  which_particle_update = dnest_get_which_particle_update();

  // only update continuum reconstruction when the corresponding parameters are updated
  if(which_parameter_update >= num_params_blr)
  {
    con_q = con_q_particles_perturb[which_particle_update];
    Fcon = Fcon_particles_perturb[which_particle_update];
    //calculate_con_from_model(model + num_params_blr*sizeof(double));
    calculate_con_from_model_semiseparable(model + num_params_blr*sizeof(double));
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  }
  else /* continuum has no change, use the previous values */
  {
    con_q = con_q_particles[which_particle_update];
    Fcon = Fcon_particles[which_particle_update];
    gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  }
  
   // only update transfer function when BLR model is changed.
   // or when forced to updated
  if( (which_parameter_update < num_params_blr-1) || force_update == 1)
  {
    TransTau = TransTau_particles_perturb[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles_perturb[which_particle_update];
    transfun_2d_cal(model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
  }
  else
  {
    TransTau = TransTau_particles[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
  }

  /* no need to calculate line when only systematic error parameter of line are updated.
   * otherwise, always need to calculate line.
   */
  if( which_parameter_update != num_params_blr-1 || force_update == 1 )
  {
    Fline2d_at_data = Fline_at_data_particles_perturb[which_particle_update];
    calculate_line2d_from_blrmodel(model, Tline_data, Vline_data_ext, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data_ext);

    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean*line_error_mean;
    for(i=0; i<n_line_data; i++)
    {
      for(j=0; j<n_vel_data; j++)
      {
        //note mask with error < 0.0
        if(Flerrs2d_data[i] > 0.0)
        {
          dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
          var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
          prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
  }
  else
  {
    /* re-point */
    Fline2d_at_data = Fline_at_data_particles[which_particle_update];
    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean*line_error_mean;
    for(i=0; i<n_line_data; i++)
    {
      for(j=0; j<n_vel_data; j++)
      {
        //note mask with error < 0.0
        if(Flerrs2d_data[i] > 0.0)
        {
          dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
          var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
          prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
  }

  return prob_line;
}
