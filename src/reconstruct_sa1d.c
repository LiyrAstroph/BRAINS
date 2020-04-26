/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_sa1d.c
 *  \brief reconstruct sa and 1d RM data.
 */

#ifdef SA

#include "brains.h"

void *best_model_sa1d;      /*!< best model */
void *best_model_std_sa1d;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess_sa1d()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  int num_ps, i, j, k, nc;
  double *pm, *pmstd;
  double *lag;
  void *posterior_sample, *post_model;
  double mean_lag, mean_lag_std, sum1, sum2;
  int size_of_modeltype = num_params * sizeof(double);
  
  best_model_sa1d = malloc(size_of_modeltype);
  best_model_std_sa1d = malloc(size_of_modeltype);

  if(thistask == roottask)
  {
    // initialize smoothing workspace
    char fname[200];
    FILE *fp, *fline, *fsa, *fsaline, *ftran, *fcon;

    sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);

    // get number of lines in posterior sample file
    get_posterior_sample_file(dnest_options_file, posterior_sample_file);

    //file for posterior sample
    fp = fopen(posterior_sample_file, "r");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file);
      exit(0);
    }
    
    fcon = fopen("data/con_rec.txt", "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/con_rec.txt.\n");
      exit(0);
    }
    //file for line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/line_rec.txt");
    fline = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }
    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran_rec.txt");
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
    nc = 0;
    
    Fcon_rm = Fcon_rm_particles[which_particle_update];
    TransTau = TransTau_particles[which_particle_update];
    Trans1D = Trans1D_particles[which_particle_update];

    Fline_sa = Fline_sa_particles[which_particle_update];
    phase_sa = phase_sa_particles[which_particle_update];
    
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
      
      calculate_con_from_model_semiseparable(post_model + num_params_blr_tot *sizeof(double));
      calculate_con_rm(post_model);
      gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);

      calculate_sa_transfun_from_blrmodel(post_model, 0);
      calculate_line_from_blrmodel(post_model, Tline, Fline, parset.n_line_recon);
      
      //if( i % (num_ps/10+1) == 0)  
      {
        for(j=0; j<parset.n_con_recon; j++)
        {
          fprintf(fcon, "%e %e %e\n", Tcon[j]*(1.0+parset.redshift), Fcon[j]/con_scale, Fcerrs[i]/con_scale);
        }
        fprintf(fcon, "\n");

        for(j=0; j<parset.n_line_recon; j++)
        {
          fprintf(fline, "%e %e\n", Tline[j]*(1.0+parset.redshift), Fline[j]/line_scale);
        }
        fprintf(fline, "\n");

        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran, "%e %e\n", TransTau[j], Trans1D[j]);
        }
        fprintf(ftran, "\n");

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
            fprintf(fsa, "%e %e\n", wave_sa_data[j], phase_sa[k*n_vel_sa_data + j]/(PhaseFactor * wave_sa_data[j]) );
          }
        }
        fprintf(fsa, "\n");
      }
    }

    fclose(fp);
    fclose(fsaline);
    fclose(fsa);
    fclose(fcon);
    fclose(ftran);
    fclose(fline);
    sa_smooth_end();

    pm = (double *)best_model_sa1d;
    pmstd = (double *)best_model_std_sa1d;
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
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_sa1d + j), 
                                             *((double *)best_model_std_sa1d+j) ); 
 
    free(post_model);
    free(posterior_sample);
  } 
}

void reconstruct_sa1d()
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
  strcat(argv[argc++], "/data/restartsa1d_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartsa1d_dnest.txt");
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
  
  reconstruct_sa1d_init();
  
  sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);
  dnest_sa1d(argc, argv);
  sa_smooth_end();

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_sa1d();

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
      Trans1D = Trans1D_particles[which_particle_update];
      
      Fline_sa = Fline_sa_particles[which_particle_update];
      phase_sa = phase_sa_particles[which_particle_update];

      sa_smooth_init(n_vel_sa_data, vel_sa_data, parset.sa_InstRes);

      //calculate_con_from_model(best_model_line1d + num_params_blr *sizeof(double));
      calculate_con_from_model_semiseparable(best_model_sa1d + num_params_blr_tot *sizeof(double));
      calculate_con_rm(best_model_sa1d);
      gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);

      calculate_sa_transfun_from_blrmodel(best_model_sa1d, 1);
      calculate_line_from_blrmodel(best_model_sa1d, Tline, Fline, parset.n_line_recon);

      sa_smooth_end();
      
      // output continuum light curve
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

      // output reconstructed line light curve
      sprintf(fname, "%s/%s", parset.file_dir, parset.pline_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<parset.n_line_recon; i++)
      {
        fprintf(fp, "%e %e\n", Tline[i]*(1.0+parset.redshift), Fline[i] / line_scale);
      }
      fclose(fp);

      // output transfer function.
      sprintf(fname, "%s/%s", parset.file_dir, parset.tran_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(i=0; i<parset.n_tau; i++)
      {
        fprintf(fp, "%e %e\n", TransTau[i], Trans1D[i]);
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
          fprintf(fp, "%e %e\n", wave_sa_data[j], phase_sa[k*n_vel_sa_data + j]/(PhaseFactor * wave_sa_data[j]) );
        }
        fprintf(fp, "\n");
      }
      fclose(fp);
    }
  }

  reconstruct_sa1d_end();

  //clear up argv
  for(i=0; i<9; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sa1d_init()
{
  int i, j;
  double dT, Tspan;

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "src/OPTIONSSA1D");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  /* RM setup */
  Tspan = (Tcon_data[n_con_data-1] - Tcon_data[0]);

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

  // set time grid for line
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Flerrs = malloc(parset.n_line_recon * sizeof(double));

  Tline_min = Tline_data[0] - fmin(0.05*(Tline_data[n_line_data - 1] - Tline_data[0]), 10.0);
  if(parset.time_back <= 0.0)
    Tline_min = fmax(Tline_min, Tcon_min + time_back_set);
  

  Tline_max = Tline_data[n_line_data -1] + fmin(0.05*(Tline_data[n_line_data - 1] - Tline_data[0]), 10.0);
  Tline_max = fmin(Tline_max, Tcon_max - 1.0);  /* The time span should be smaller than that of the continuum */
  
  if(thistask == roottask)
    printf("Tline_min_max: %f %f\n", Tline_min-Tline_data[0], Tline_max - Tline_data[n_line_data-1]);

  dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

  for(i=0; i<parset.n_line_recon; i++)
  {
    Tline[i] = Tline_min + i*dT;
  }

  // Fcon perturbed and accepted for each particle
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

  Trans1D_particles = malloc(parset.num_particles * sizeof(double *));
  Trans1D_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans1D_particles[i] = malloc(parset.n_tau * sizeof(double));
    Trans1D_particles_perturb[i] = malloc(parset.n_tau * sizeof(double));
  }

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_line_data * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_line_data * sizeof(double));
  }

  /* SA setup */
  phase_sa_particles = malloc(parset.num_particles * sizeof(double *));
  phase_sa_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    phase_sa_particles[i] = malloc(n_base_sa_data * n_vel_sa_data * sizeof(double));
    phase_sa_particles_perturb[i] = malloc(n_base_sa_data * n_vel_sa_data * sizeof(double));
  }

  Fline_sa_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_sa_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_sa_particles[i] = malloc(n_vel_sa_data * sizeof(double));
    Fline_sa_particles_perturb[i] = malloc(n_vel_sa_data * sizeof(double));
  }

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
    if(parset.n_cloud_per_task <= 1000)
      icr_cloud_save = 1;
    else
      icr_cloud_save = parset.n_cloud_per_task/1000;

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

void reconstruct_sa1d_end()
{
  int i;
  
  /* RM setup */
  free(Tline);
  free(Fline);
  free(Flerrs);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Fcon_rm_particles[i]);
    free(Fcon_rm_particles_perturb[i]);
  }
  free(Fcon);
  free(Fcon_rm_particles);
  free(Fcon_rm_particles_perturb);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Trans1D_particles[i]);
    free(Trans1D_particles_perturb[i]);
    free(TransTau_particles[i]);
    free(TransTau_particles_perturb[i]);
    free(Fline_at_data_particles[i]);
    free(Fline_at_data_particles_perturb[i]);
  }
  free(Trans1D_particles);
  free(Trans1D_particles_perturb);
  free(TransTau_particles);
  free(TransTau_particles_perturb);
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);

  /* SA setup */
  for(i=0; i<parset.num_particles; i++)
  {
    free(phase_sa_particles[i]);
    free(phase_sa_particles_perturb[i]);

    free(Fline_sa_particles[i]);
    free(Fline_sa_particles_perturb[i]);
  }
  free(phase_sa_particles);
  free(phase_sa_particles_perturb);

  free(Fline_sa_particles);
  free(Fline_sa_particles_perturb);

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
  free(best_model_sa1d);
  free(best_model_std_sa1d);

  free(clouds_weight);
  free(clouds_alpha);
  free(clouds_beta);
  free(clouds_tau);

  free(workspace_phase);

  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }
  if(thistask == roottask)
  {
    printf("Ends reconstruct_sa1d.\n");
  }
  return;
}

/*!
 * this function calculate probability.
 * 
 * At each MCMC step, only one parameter is updated, which only changes some values; thus,
 * optimization that reuses the unchanged values can improve computation efficiency.
 */
double prob_sa1d(const void *model)
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

  /* only update transfer function when BLR model is changed
   * or forced to update (force_update = 1)
   * Trans1D is a pointer to the transfer function
   */
  if(     ((which_parameter_update < num_params_blr_tot)
        && (which_parameter_update != num_params_blr-1)) 
        || force_update == 1)
  {
    /* re-point */
    TransTau = TransTau_particles_perturb[which_particle_update];
    Trans1D = Trans1D_particles_perturb[which_particle_update]; 
    
    phase_sa = phase_sa_particles_perturb[which_particle_update];
    Fline_sa = Fline_sa_particles_perturb[which_particle_update];

    calculate_sa_transfun_from_blrmodel(model, 0);

    /* caclulate prob_sa */
    for(j=0; j<n_epoch_sa_data; j++)
    {
      for(i=0; i<n_vel_sa_data; i++)
      {
        dy = Fline_sa[i] - Fline_sa_data[i + j*n_vel_sa_data];
        var2 = Flerrs_sa_data[i+j*n_vel_sa_data]*Flerrs_sa_data[i+j*n_vel_sa_data];
        prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
      }
    }
    for(j=0; j<n_base_sa_data; j++)
    {
      for(i=0; i<n_vel_sa_data; i++)
      {
        dy = phase_sa[j*n_vel_sa_data + i] - phase_sa_data[j*n_vel_sa_data + i];
        var2 = pherrs_sa_data[j*n_vel_sa_data + i] * pherrs_sa_data[j*n_vel_sa_data + i];
        prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
      }
    }
    prob_sa_particles_perturb[which_particle_update] = prob_sa;
  }
  else
  {
    /* re-point */
    TransTau = TransTau_particles[which_particle_update];
    Trans1D = Trans1D_particles[which_particle_update];

    phase_sa = phase_sa_particles[which_particle_update];
    Fline_sa = Fline_sa_particles[which_particle_update];

    prob_sa = prob_sa_particles[which_particle_update];
  }

  /* no need to calculate line when only systematic error parameter of line are updated.
   * otherwise, always need to calculate line.
   */
  if( which_parameter_update != num_params_blr-1 || force_update == 1 )
  {
    /* re-point */
    Fline_at_data = Fline_at_data_particles_perturb[which_particle_update];
    calculate_line_from_blrmodel(model, Tline_data, Fline_at_data, n_line_data);

    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
    for(i=0; i<n_line_data; i++)
    {
      //note mask with error < 0.0
      if(Flerrs_data[i] > 0.0)
      {
        dy = Fline_data[i] - Fline_at_data[i] ;
        var2 = Flerrs_data[i]*Flerrs_data[i] + var2_se;
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }
  else
  {
    /* re-point */
    Fline_at_data = Fline_at_data_particles[which_particle_update];
    var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
    for(i=0; i<n_line_data; i++)
    {
      //note mask with error < 0.0
      if(Flerrs_data[i] > 0.0)
      {
        dy = Fline_data[i] - Fline_at_data[i] ;
        var2 = Flerrs_data[i]*Flerrs_data[i] + var2_se;
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
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
double prob_initial_sa1d(const void *model)
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
  
  TransTau = TransTau_particles[which_particle_update];
  Trans1D = Trans1D_particles[which_particle_update]; 
    
  phase_sa = phase_sa_particles[which_particle_update];
  Fline_sa = Fline_sa_particles[which_particle_update];

  calculate_sa_transfun_from_blrmodel(model, 0);

  /* caclulate prob_sa */
  for(j=0; j<n_epoch_sa_data; j++)
  {
    for(i=0; i<n_vel_sa_data; i++)
    {
      dy = Fline_sa[i] - Fline_sa_data[i + j*n_vel_sa_data];
      var2 = Flerrs_sa_data[i+j*n_vel_sa_data]*Flerrs_sa_data[i+j*n_vel_sa_data];
      prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
    }
  }
  for(j=0; j<n_base_sa_data; j++)
  {
    for(i=0; i<n_vel_sa_data; i++)
    {
      dy = phase_sa[j*n_vel_sa_data + i] - phase_sa_data[j*n_vel_sa_data + i];
      var2 = pherrs_sa_data[j*n_vel_sa_data + i] * pherrs_sa_data[j*n_vel_sa_data + i];
      prob_sa += -0.5 * (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI);
    }
  }
  prob_sa_particles[which_particle_update] = prob_sa;
  
  Fline_at_data = Fline_at_data_particles_perturb[which_particle_update];
  calculate_line_from_blrmodel(model, Tline_data, Fline_at_data, n_line_data);
  var2_se = (exp(pm[num_params_blr-1])-1.0) * (exp(pm[num_params_blr-1])-1.0) * line_error_mean_sq;
  for(i=0; i<n_line_data; i++)
  {
    //note mask with error < 0.0
    if(Flerrs_data[i] > 0.0)
    {
      dy = Fline_data[i] - Fline_at_data[i] ;
      var2 = Flerrs_data[i]*Flerrs_data[i] + var2_se;
      prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
    }
  }
  
  return prob_sa + prob_line;
}
#endif