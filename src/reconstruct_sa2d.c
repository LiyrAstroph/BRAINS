/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_sa2d.c
 *  \brief reconstruct sa and 2d RM data.
 */

#ifdef SpecAstro

#include "brains.h"

void *best_model_sa2d;      /*!< best model */
void *best_model_std_sa2d;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess_sa2d()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  int num_ps, i, j, k;
  double *pm, *pmstd;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);
  
  best_model_sa2d = malloc(size_of_modeltype);
  best_model_std_sa2d = malloc(size_of_modeltype);

  if(thistask == roottask)
  {
    // initialize smoothing workspace
    smooth_init(n_vel_data_ext, Vline_data_ext);
    sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);

    char fname[200];
    FILE *fp, *fline, *fline1d, *fsa, *fsaline, *ftran, *fcon, *fcon_rm;
    double *Fline1d, dV;
    
    // velocity grid width, in term of wavelength of Hbeta.
    dV = (Vline_data[n_vel_data-1]-Vline_data[0])/(n_vel_data-1) * parset.linecenter/C_Unit; 
    Fline1d = malloc(n_line_data * sizeof(double));

    // get number of lines in posterior sample file
    //get_posterior_sample_file(dnest_options_file, posterior_sample_file);
    dnest_get_posterior_sample_file(posterior_sample_file);

    //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    sprintf(fname, "%s/%s", parset.file_dir, "data/con_rec.txt");
    fcon = fopen(fname, "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/con_rec.txt.\n");
      exit(0);
    }

    //file for continuum reconstruction after detrending
    if(parset.flag_trend_diff > 0)
    {
      sprintf(fname, "%s/%s", parset.file_dir, "data/con_rm_rec.txt");
      fcon_rm = fopen(fname, "w");
      if(fcon_rm == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file data/con_rm_rec.txt.\n");
        exit(0);
      }
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

    //file for sa line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/sa_line_rec.txt");
    fsaline = fopen(fname, "w");
    if(fsaline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for phase
    sprintf(fname, "%s/%s", parset.file_dir, "data/sa_phase_rec.txt");
    fsa = fopen(fname, "w");
    if(fsa == NULL)
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

    post_model = malloc(size_of_modeltype);
    posterior_sample = malloc(num_ps * size_of_modeltype);

    force_update = 1;
    which_parameter_update = -1; // force to update the transfer function
    which_particle_update = 0;

    Fcon_rm = Fcon_rm_particles[which_particle_update];
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
      calculate_con_from_model_semiseparable(post_model + num_params_blr_tot *sizeof(double));
      calculate_con_rm(post_model);
      gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);

      calculate_sa_transfun2d_from_blrmodel(post_model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
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

      //if( i % (num_ps/10+1) == 0)  
      {
        // output continuum
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%e %e\n", Tcon[j]*(1.0+parset.redshift), Fcon[j]/con_scale);
        }
        fprintf(fcon, "\n");

        if(parset.flag_trend_diff > 0)
        {
          for(j=0; j<parset.n_con_recon; j++)
          {
            fprintf(fcon_rm, "%e %e\n", Tcon[j]*(1.0+parset.redshift), Fcon_rm[j]/con_scale);
          }
          fprintf(fcon_rm, "\n");
        }      

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

        // output sa line
        for(j=0; j<n_vel_sa_data; j++)
        {
          fprintf(fsaline, "%e %e\n", wave_sa_data[j], Fline_sa[j]);
        }
        fprintf(fsaline, "\n");

        // output sa phase
        for(k=0; k<n_base_sa_data; k++)
        {
          for(j=0; j<n_vel_sa_data; j++)
          {
            fprintf(fsa, "%e %e\n", wave_sa_data[j], phase_sa[k*n_vel_sa_data + j]/(ScaleFactor[j]) );
          }
        }
        fprintf(fsa, "\n");
      }
    }

    smooth_end();
    sa_smooth_end();
    fclose(fp);
    fclose(fcon);
    if(parset.flag_trend_diff > 0)
      fclose(fcon_rm);
    fclose(fline);
    fclose(fline1d);
    fclose(ftran);
    fclose(fsaline);
    fclose(fsa);
    
    pm = (double *)best_model_sa2d;
    pmstd = (double *)best_model_std_sa2d;
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

    for(j = 0; j<num_params_blr_tot + num_params_var; j++)
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_sa2d + j), 
                                             *((double *)best_model_std_sa2d+j) ); 

    free(post_model);
    free(posterior_sample);
    free(Fline1d);
  }
}

void reconstruct_sa2d()
{
  int i, argc=0;
  char **argv;

  //configure restart of dnest
  argv = malloc(11*sizeof(char *));
  for(i=0; i<11; i++)
  {
    argv[i] = malloc(BRAINS_MAX_STR_LENGTH*sizeof(char));
  }
  //setup argc and argv
  strcpy(argv[argc++], "dnest");
  strcpy(argv[argc++], "-s");
  strcpy(argv[argc], parset.file_dir);
  strcat(argv[argc++], "/data/restartsa2d_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartsa2d_dnest.txt");
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
  
  // sample tag 
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "_sa2d");

  reconstruct_sa2d_init();
  
  smooth_init(n_vel_data_ext, Vline_data_ext);
  sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);
  dnest_sa2d(argc, argv);
  smooth_end();
  sa_smooth_end();
 
  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_sa2d();
    if(thistask == roottask)
    {
      FILE *fp;
      char fname[200];
      int j, k;

      force_update = 1;
      which_parameter_update = -1; // force to update the transfer function
      which_particle_update = 0;

      Fcon_rm = Fcon_rm_particles[which_particle_update];
      TransTau = TransTau_particles[which_particle_update];
      Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
      Fline2d_at_data = Fline_at_data_particles[which_particle_update];

      //calculate_con_from_model(best_model_line2d + num_params_blr *sizeof(double));
      calculate_con_from_model_semiseparable(best_model_sa2d + num_params_blr_tot *sizeof(double));
      calculate_con_rm(best_model_sa2d);
      gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);

      sprintf(fname, "%s/%s", parset.file_dir, parset.pcon_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      for(i=0; i<parset.n_con_recon; i++)
      {
        fprintf(fp, "%e %e %e\n", Tcon[i]*(1.0+parset.redshift), Fcon[i] / con_scale, Fcerrs[i] / con_scale);
      }
      fclose(fp);
      
      smooth_init(n_vel_data_ext, Vline_data_ext);
      sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);

      calculate_sa_transfun2d_from_blrmodel(best_model_sa2d, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
      calculate_line2d_from_blrmodel(best_model_sa2d, Tline_data, Vline_data_ext, Trans2D_at_veldata, 
                                      Fline2d_at_data, n_line_data, n_vel_data_ext);
      
      sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_data_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
    
      fprintf(fp, "# %d %d\n", n_line_data, n_vel_data);
      for(i=0; i<n_line_data; i++)
      {
        fprintf(fp, "# %f\n", Tline_data[i]*(1.0+parset.redshift));
        for(j=0; j<n_vel_data; j++)
        {
          fprintf(fp, "%e %e\n", Wline_data[j],   
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
        fprintf(fp, "# %f\n", TransTau[i]);
        for(j=0; j<n_vel_data; j++)
        {
          fprintf(fp, "%e %e\n", Vline_data[j]*VelUnit, Trans2D_at_veldata[i*n_vel_data_ext + (j+n_vel_data_incr)]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      // output sa line 
      sprintf(fname, "%s/%s", parset.file_dir, "data/psa_line.txt");
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
        exit(0);
      }
      for(j=0; j<n_vel_sa_data; j++)
      {
        fprintf(fp, "%e %e\n", wave_sa_data[j], Fline_sa[j]);
      }
      fclose(fp);

      //output sa phase
      sprintf(fname, "%s/%s", parset.file_dir, "data/psa_phase.txt");
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
        exit(0);
      }
      for(k=0; k<n_base_sa_data; k++)
      {
        for(j=0; j<n_vel_sa_data; j++)
        {
          fprintf(fp, "%e %e\n", wave_sa_data[j], phase_sa[k*n_vel_sa_data + j]/(ScaleFactor[j]) );
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      smooth_end();

      // recovered line2d at specified points
      smooth_init(parset.n_vel_recon, TransV);

      which_parameter_update = -1;
      which_particle_update = 0;

      transfun_2d_cal(best_model_sa2d, TransV, Trans2D, parset.n_vel_recon, 0);

      /* there is no data for spectral broadening at given specified epoch, using the mean value
       * and set InstRes_err=0.0.
       */ 
      double *pm = (double *)best_model_sa2d;
      if(parset.flag_InstRes > 1) 
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
          linecenter_mean += pm[idx_linecenter + i];
        }
        pm[idx_linecenter] = linecenter_mean/n_line_data;
      }
      
      calculate_line2d_from_blrmodel(best_model_sa2d, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

      sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }

      fprintf(fp, "# %d %d\n", parset.n_line_recon, parset.n_vel_recon);
      for(i=0; i<parset.n_line_recon; i++)
      {
        fprintf(fp, "# %f\n", Tline[i]*(1.0+parset.redshift));
        for(j=0; j<parset.n_vel_recon; j++)
        {
          //fprintf(fp, "%e %e\n", TransV[j]*VelUnit,  Fline2d[i*parset.n_vel_recon + j] / line_scale);
          fprintf(fp, "%e %e\n", (1.0 + TransV[j]/C_Unit) * parset.linecenter * (1.0+parset.redshift),  
                Fline2d[i*parset.n_vel_recon + j] / line_scale);
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
        fprintf(fp, "# %f\n", TransTau[i]);
        for(j=0; j<parset.n_vel_recon; j++)
        {
          fprintf(fp, "%e %e\n", TransV[j]*VelUnit, Trans2D[i*parset.n_vel_recon + j]);
        }

        fprintf(fp, "\n");
      }
      fclose(fp);

      smooth_end();
      sa_smooth_end();
    }
  }
  
  reconstruct_sa2d_end();

  //clear up argv
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sa2d_init()
{
  int i, j;
  double dT, Tspan;

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "param/OPTIONSSA2D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);
  
  /* RM setup */
  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  Tcon_min = Tcon_data[0] - time_back_set - 10.0;
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 20.0);
  Tcon_max = fmax(Tcon_max, Tline_data[n_line_data -1] + 10.0);  /* The time span should cover that of the emission line data */
  
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

  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));
  
  Tline_min = Tline_data[0] - fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10.0);
  if(parset.time_back <= 0.0)
    Tline_min = fmax(Tline_min, Tcon_min + time_back_set);

  Tline_max = Tline_data[n_line_data -1] + fmin(0.1*(Tline_data[n_line_data - 1] - Tline_data[0]), 10.0);
  Tline_max = fmin(Tline_max, Tcon_max - 1.0);  /* The time span should be smaller than that of the continuum */

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

  Fcon = malloc(parset.n_con_recon * sizeof(double *));
  Fcon_rm_particles = malloc(parset.num_particles * sizeof(double *));
  Fcon_rm_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fcon_rm_particles[i] = malloc(parset.n_con_recon * sizeof(double));
    Fcon_rm_particles_perturb[i] = malloc(parset.n_con_recon * sizeof(double));
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

  /* SA setup */
  phase_sa = malloc(n_base_sa_data * n_vel_sa_data * sizeof(double));
  Fline_sa = malloc(n_vel_sa_data * sizeof(double));

  prob_sa_particles = malloc(parset.num_particles * sizeof(double));
  prob_sa_particles_perturb = malloc(parset.num_particles * sizeof(double));

  /* cloud sample related */
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));

  workspace_phase = malloc( (3*n_vel_sa_data)* sizeof(double));

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
}

void reconstruct_sa2d_end()
{
  int i;

  /* RM setup */
  free(Tline);
  free(Fline2d);
  free(TransV);
  free(Trans2D);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Fcon_rm_particles[i]);
    free(Fcon_rm_particles_perturb[i]);
  }
  free(Fcon);
  free(Fcon_rm_particles);
  free(Fcon_rm_particles_perturb);
  
  free(Fline_at_data);
  
  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans2D_at_veldata_particles[i]);
    free(Trans2D_at_veldata_particles_perturb[i]);
    free(TransTau_particles[i]);
    free(TransTau_particles_perturb[i]);
    free(Fline_at_data_particles[i]);
    free(Fline_at_data_particles_perturb[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(Trans2D_at_veldata_particles_perturb);
  free(TransTau_particles);
  free(TransTau_particles_perturb);
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);

  /* SA setup */
  free(phase_sa);
  free(Fline_sa);

  free(prob_sa_particles);
  free(prob_sa_particles_perturb);

  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
    free(par_prior_gaussian[i]);
  }
  free(par_range_model);
  free(par_prior_gaussian);
  free(par_prior_model);

  free(par_fix);
  free(par_fix_val);
  free(best_model_sa2d);
  free(best_model_std_sa2d);

  free(clouds_weight);
  free(clouds_alpha);
  free(clouds_beta);
  free(clouds_tau);
  free(clouds_vel);

  free(workspace_phase);

  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  if(thistask == roottask)
  {
    printf("Ends reconstruct_sa2d.\n");
  }
}

/*!
 * this function calculate probability.
 * 
 * At each MCMC step, only one parameter is updated, which only changes some values; thus,
 * optimization that reuses the unchanged values can improve computation efficiency.
 */
double prob_sa2d(const void *model)
{
  double prob_sa = 0.0, prob_line=0.0, var2, dy, var2_se;
  int i, j;
  double *pm = (double *)model;
  
  which_particle_update = dnest_get_which_particle_update();

  // only update continuum reconstruction when the corresponding parameters are updated
  if( which_parameter_update >= num_params_blr_tot )
  {
    Fcon_rm = Fcon_rm_particles_perturb[which_particle_update];
    //calculate_con_from_model(model + num_params_blr*sizeof(double));
    calculate_con_from_model_semiseparable(model + num_params_blr_tot*sizeof(double));
    calculate_con_rm(model);

    gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);
  }
  else /* continuum has no change, use the previous values */
  {
    Fcon_rm = Fcon_rm_particles[which_particle_update];
    gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);
  }

  /* only update transfer function when BLR model of RM and SA is changed
   * or forced to update (force_update = 1)
   * Trans1D is a pointer to the transfer function
   */
  if((
      (which_parameter_update < num_params_blr_model) || 
      (which_parameter_update >= num_params_blr && which_parameter_update < num_params_blr_tot)
     ) || force_update == 1 )
  {
    /* re-point */
    TransTau = TransTau_particles_perturb[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles_perturb[which_particle_update];

    calculate_sa_transfun2d_from_blrmodel(model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
    
    /* caclulate prob_sa */
    for(j=0; j<n_epoch_sa_data; j++)
    {
      for(i=0; i<n_vel_sa_data; i++)
      {
        //note mask with error < 0.0
        if(Flerrs_sa_data[i+j*n_vel_sa_data] > 0.0)
        {
          dy = Fline_sa[i] - Fline_sa_data[i + j*n_vel_sa_data];
          var2 = Flerrs_sa_data[i+j*n_vel_sa_data]*Flerrs_sa_data[i+j*n_vel_sa_data];
          prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
        }
      }
    }

    for(j=0; j<n_base_sa_data; j++)
    {
      for(i=0; i<n_vel_sa_data; i++)
      {
        //note mask with error < 0.0
        if(pherrs_sa_data[j*n_vel_sa_data + i] > 0.0)
        {
          dy = phase_sa[j*n_vel_sa_data + i] - phase_sa_data[j*n_vel_sa_data + i];
          var2 = pherrs_sa_data[j*n_vel_sa_data + i] * pherrs_sa_data[j*n_vel_sa_data + i];
          prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
    prob_sa_particles_perturb[which_particle_update] = prob_sa;
  }
  else
  {
    /* re-point */
    TransTau = TransTau_particles[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];

    prob_sa = prob_sa_particles[which_particle_update];
  }
  
  /* no need to calculate line when only systematic error parameter of line are updated.
   * otherwise, always need to calculate line.
   */
  if( which_parameter_update != num_params_blr-1 || force_update == 1 )
  {
    /* re-point */
    Fline2d_at_data = Fline_at_data_particles_perturb[which_particle_update];
    calculate_line2d_from_blrmodel(model, Tline_data, Vline_data_ext, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data_ext);
    
    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
    for(i=0; i<n_line_data; i++)
    {
      for(j=0; j<n_vel_data; j++)
      {
        //note mask with error < 0.0
        if(Flerrs2d_data[i*n_vel_data+j] > 0.0)
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
    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
    for(i=0; i<n_line_data; i++)
    {
      for(j=0; j<n_vel_data; j++)
      {
        //note mask with error < 0.0
        if(Flerrs2d_data[i*n_vel_data+j] > 0.0)
        {
          dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
          var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
          prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
  }

  return prob_sa + prob_line;
}

/*!
 * this function calculate probability.
 * 
 * At each MCMC step, only one parameter is updated, which only changes some values; thus,
 * optimization that reuses the unchanged values can improve computation efficiency.
 */
double prob_initial_sa2d(const void *model)
{
  double prob_sa = 0.0, prob_line=0.0, var2, dy, var2_se;
  int i, j;
  double *pm = (double *)model;
  
  which_particle_update = dnest_get_which_particle_update();

  Fcon_rm = Fcon_rm_particles[which_particle_update];
  //calculate_con_from_model(model + num_params_blr*sizeof(double));
  calculate_con_from_model_semiseparable(model + num_params_blr_tot*sizeof(double));
  calculate_con_rm(model);
  gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);
  

  /* re-point */
  TransTau = TransTau_particles[which_particle_update];
  Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];

  calculate_sa_transfun2d_from_blrmodel(model, Vline_data_ext, Trans2D_at_veldata, n_vel_data_ext, 0);
    
  /* caclulate prob_sa */
  for(j=0; j<n_epoch_sa_data; j++)
  {
    for(i=0; i<n_vel_sa_data; i++)
    {
      //note mask with error < 0.0
      if(Flerrs_sa_data[i+j*n_vel_sa_data] > 0.0)
      {
        dy = Fline_sa[i] - Fline_sa_data[i + j*n_vel_sa_data];
        var2 = Flerrs_sa_data[i+j*n_vel_sa_data]*Flerrs_sa_data[i+j*n_vel_sa_data];
        prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
      }
    }
  }
  for(j=0; j<n_base_sa_data; j++)
  {
    for(i=0; i<n_vel_sa_data; i++)
    {
      //note mask with error < 0.0
      if(pherrs_sa_data[j*n_vel_sa_data + i] > 0.0)
      {
        dy = phase_sa[j*n_vel_sa_data + i] - phase_sa_data[j*n_vel_sa_data + i];
        var2 = pherrs_sa_data[j*n_vel_sa_data + i] * pherrs_sa_data[j*n_vel_sa_data + i];
        prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
      }
    }
  }
  prob_sa_particles[which_particle_update] = prob_sa;
  
  Fline2d_at_data = Fline_at_data_particles[which_particle_update];
  calculate_line2d_from_blrmodel(model, Tline_data, Vline_data_ext, Trans2D_at_veldata, Fline2d_at_data, n_line_data, n_vel_data_ext);
  
  var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
  for(i=0; i<n_line_data; i++)
  {
    for(j=0; j<n_vel_data; j++)
    {
      //note mask with error < 0.0
      if(Flerrs2d_data[i*n_vel_data+j] > 0.0)
      {
        dy = Fline2d_data[i*n_vel_data + j] - Fline2d_at_data[i * n_vel_data_ext + (j+n_vel_data_incr)];
        var2 = Flerrs2d_data[i*n_vel_data+j]*Flerrs2d_data[i*n_vel_data+j] + var2_se;
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }

  return prob_sa + prob_line;
}

#endif