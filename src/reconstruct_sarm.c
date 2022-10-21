/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */


#ifdef SpecAstro

#include "brains.h"

void *best_model_sarm;      /*!< best model */
void *best_model_std_sarm;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocess_sarm()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  int num_ps, i, j, k, m;
  double *pm, *pmstd;
  void *posterior_sample, *post_model;
  int size_of_modeltype = num_params * sizeof(double);

  best_model_sarm = malloc(size_of_modeltype);
  best_model_std_sarm = malloc(size_of_modeltype);

  if(thistask == roottask)
  {
    char fname[200];
    FILE *fp, *fcon, *fcon_rm, *fline, *ftran, *ftran_alpha, *ftran_beta, *fsa;

    sarm_smooth_init(n_vel_sarm_data_ext, vel_sa_data_ext, parset.sa_InstRes);
    
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
    fcon = fopen("data/con_rec.txt", "w");
    if(fcon == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file data/con_rec.txt.\n");
      exit(0);
    }

    //file for continuum reconstruction after detrending
    if(parset.flag_trend_diff > 0)
    {
      fcon_rm = fopen("data/con_rm_rec.txt", "w");
      if(fcon_rm == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file data/con_rm_rec.txt.\n");
        exit(0);
      }
    }
    
    //file for line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/line2d_sarm_rec.txt");
    fline = fopen(fname, "w");
    if(fline == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran2d_sarm_rec.txt");
    ftran = fopen(fname, "w");
    if(ftran == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran2d_alpha_sarm_rec.txt");
    ftran_alpha = fopen(fname, "w");
    if(ftran == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for transfer function
    sprintf(fname, "%s/%s", parset.file_dir, "data/tran2d_beta_sarm_rec.txt");
    ftran_beta = fopen(fname, "w");
    if(ftran == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
      exit(0);
    }

    //file for phase
    sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_phase_rec.txt");
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
    Trans_alpha_at_veldata = Trans_alpha_at_veldata_particles[which_particle_update];
    Trans_beta_at_veldata = Trans_beta_at_veldata_particles[which_particle_update];
    Fline2d_sarm_at_data = Fline_at_data_particles[which_particle_update];
    phase_sarm_at_data = phase_at_data_particles[which_particle_update];

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

      transfun_sarm_cal_cloud(post_model, vel_sa_data_ext, Trans2D_at_veldata, Trans_alpha_at_veldata,
                            Trans_beta_at_veldata, n_vel_sarm_data_ext, 0);
  
      calculate_sarm_with_sample(post_model);

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
        for(j=0; j<n_epoch_sarm_data; j++)
        {
          for(k=0; k<n_vel_sarm_data; k++)
          {
            fprintf(fline, "%e ", Fline2d_sarm_at_data[j * n_vel_sarm_data_ext + (k+n_vel_sa_data_incr)]/line_sarm_scale);
          }
          fprintf(fline, "\n");
        }
        fprintf(fline, "\n");

        // output transfer function
        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran, "%e ", TransTau[j]);
          
          for(k=0; k<n_vel_sarm_data; k++)
          {
            fprintf(ftran, "%e ", Trans2D_at_veldata[j * n_vel_sarm_data_ext + (k+n_vel_sa_data_incr)]);
          }
          fprintf(ftran, "\n");
        }
        fprintf(ftran, "\n"); 

        // output transfer function alpha
        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran_alpha, "%e ", TransTau[j]);
          
          for(k=0; k<n_vel_sarm_data; k++)
          {
            fprintf(ftran_alpha, "%e ", Trans_alpha_at_veldata[j * n_vel_sarm_data_ext + (k+n_vel_sa_data_incr)]);
          }
          fprintf(ftran_alpha, "\n");
        }
        fprintf(ftran_alpha, "\n"); 

        // output transfer function beta
        for(j=0; j<parset.n_tau; j++)
        {
          fprintf(ftran_beta, "%e ", TransTau[j]);
          
          for(k=0; k<n_vel_sarm_data; k++)
          {
            fprintf(ftran_beta, "%e ", Trans_beta_at_veldata[j * n_vel_sarm_data_ext + (k+n_vel_sa_data_incr)]);
          }
          fprintf(ftran_beta, "\n");
        }
        fprintf(ftran_beta, "\n"); 

        // output sarm phase
        for(m=0; m<n_epoch_sarm_data; m++)
        {
          for(k=0; k<n_base_sarm_data; k++)
          {
            for(j=0; j<n_vel_sarm_data; j++)
            {
              fprintf(fsa, "%e ", phase_sarm_at_data[m*n_vel_sarm_data_ext*n_base_sarm_data + k*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)]
                  /(ScaleFactor[j]) );
            }
            fprintf(fsa, "\n");
          }
          fprintf(fsa, "\n");    
        }
      }
    }

    sarm_smooth_end();

    fclose(fp);
    fclose(fcon);
    if(parset.flag_trend_diff > 0)
      fclose(fcon_rm);
    fclose(fline);
    fclose(ftran);
    fclose(ftran_alpha);
    fclose(ftran_beta);
    fclose(fsa);

    pm = (double *)best_model_sarm;
    pmstd = (double *)best_model_std_sarm;
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
      printf("Best params %d %f +- %f\n", j, *((double *)best_model_sarm + j), 
                                             *((double *)best_model_std_sarm+j) ); 

    free(post_model);
    free(posterior_sample);
  }
}

void reconstruct_sarm()
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
  strcat(argv[argc++], "/data/restartsarm_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartsarm_dnest.txt");
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
  strcpy(argv[argc++], "_sarm");

  reconstruct_sarm_init();
  sarm_smooth_init(n_vel_sarm_data_ext, vel_sa_data_ext, parset.sa_InstRes);
  
  dnest_sarm(argc, argv);

  sarm_smooth_end();

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocess_sarm();

    if(thistask == roottask)
    {
      FILE *fp;
      char fname[200];
      int j, k, m;

      force_update = 1;
      which_parameter_update = -1; // force to update the transfer function
      which_particle_update = 0;
      
      sarm_smooth_init(n_vel_sarm_data_ext, vel_sa_data_ext, parset.sa_InstRes);

      Fcon_rm = Fcon_rm_particles[which_particle_update];
      TransTau = TransTau_particles[which_particle_update];
      Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
      Trans_alpha_at_veldata = Trans_alpha_at_veldata_particles[which_particle_update];
      Trans_beta_at_veldata = Trans_beta_at_veldata_particles[which_particle_update];
      Fline2d_sarm_at_data = Fline_at_data_particles[which_particle_update];
      phase_sarm_at_data = phase_at_data_particles[which_particle_update];

      calculate_con_from_model_semiseparable(best_model_sarm + num_params_blr_tot *sizeof(double));
      calculate_con_rm(best_model_sarm);
      gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);

      transfun_sarm_cal_cloud(best_model_sarm, vel_sa_data_ext, Trans2D_at_veldata, Trans_alpha_at_veldata,
                            Trans_beta_at_veldata, n_vel_sarm_data_ext, 0);
  
      calculate_sarm_with_sample(best_model_sarm);

      sprintf(fname, "%s/%s", parset.file_dir, parset.pline2d_data_out_file);
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
    
      fprintf(fp, "# %d %d\n", n_epoch_sarm_data, n_vel_sarm_data);
      for(i=0; i<n_epoch_sarm_data; i++)
      {
        fprintf(fp, "# %f\n", Tline_sarm_data[i]*(1.0+parset.redshift));
        for(j=0; j<n_vel_sarm_data; j++)
        {
          fprintf(fp, "%e %e\n", wave_sa_data[j],   
                            Fline2d_sarm_at_data[i*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)] / line_sarm_scale);
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
    
      fprintf(fp, "# %d %d\n", parset.n_tau, n_vel_sarm_data);
      for(i=0; i<parset.n_tau; i++)
      {
        fprintf(fp, "# %f\n", TransTau[i]);
        for(j=0; j<n_vel_sarm_data; j++)
        {
          fprintf(fp, "%e %e\n", vel_sa_data[j]*VelUnit, Trans2D_at_veldata[i*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)]);
        }
        fprintf(fp, "\n");
      }
      fclose(fp);

      /* output sarm phase */
      sprintf(fname, "%s/%s", parset.file_dir, "data/psarm_phase_data.txt");
      fp = fopen(fname, "w");
      if(fp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s\n", fname);
        exit(-1);
      }
      for(m=0; m<n_epoch_sarm_data; m++)
      {
        for(k=0; k<n_base_sarm_data; k++)
        {
          for(j=0; j<n_vel_sarm_data; j++)
          {
            fprintf(fp, "%e %e\n", wave_sa_data[j], 
                phase_sarm_at_data[m*n_vel_sarm_data_ext*n_base_sarm_data + k*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)]
                /(ScaleFactor[j]) );
          }
          fprintf(fp, "\n");
        }
      }
      fclose(fp);

      sarm_smooth_end();
    }
  }

  reconstruct_sarm_end();

  //clear up argv
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_sarm_init()
{
  int i, j;
  double dT, Tspan;

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "param/OPTIONSSARM");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  /* RM setup */
  Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  Tcon_min = Tcon_data[0] - time_back_set - 10.0;
  Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 20.0);
  Tcon_max = fmax(Tcon_max, Tline_sarm_data[n_epoch_sarm_data -1] + 10.0);  /* The time span should cover that of the emission line data */

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
  double vel_max_set = vel_sa_data_ext[n_vel_sarm_data_ext -1], vel_min_set = vel_sa_data_ext[0];
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
    Trans2D_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
    Trans2D_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
  }

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_epoch_sarm_data * n_vel_sarm_data_ext * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_epoch_sarm_data * n_vel_sarm_data_ext * sizeof(double));
  }

  /* SA setup */
  Trans_alpha_at_veldata_particles = malloc(parset.num_particles * sizeof(double *));
  Trans_alpha_at_veldata_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  Trans_beta_at_veldata_particles = malloc(parset.num_particles * sizeof(double *));
  Trans_beta_at_veldata_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Trans_alpha_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
    Trans_alpha_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
    Trans_beta_at_veldata_particles[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
    Trans_beta_at_veldata_particles_perturb[i] = malloc(parset.n_tau * n_vel_sarm_data_ext * sizeof(double));
  }

  phase_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  phase_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    phase_at_data_particles[i] = malloc(n_epoch_sarm_data * n_base_sarm_data * n_vel_sarm_data_ext * sizeof(double));
    phase_at_data_particles_perturb[i] = malloc(n_epoch_sarm_data * n_base_sarm_data * n_vel_sarm_data_ext * sizeof(double));
  }

  /* cloud sample related */
  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));

  momentum_sarm_alpha = malloc(n_epoch_sarm_data * n_vel_sarm_data_ext * sizeof(double));
  momentum_sarm_beta = malloc(n_epoch_sarm_data * n_vel_sarm_data_ext * sizeof(double));
  
  workspace_phase = malloc( (3*n_vel_sarm_data_ext)* sizeof(double));

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

void reconstruct_sarm_end()
{
  int i;

  /* RM setup */
  free(TransV);

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
    
    free(Trans_alpha_at_veldata_particles[i]);
    free(Trans_alpha_at_veldata_particles_perturb[i]);
    free(Trans_beta_at_veldata_particles[i]);
    free(Trans_beta_at_veldata_particles_perturb[i]);
    free(phase_at_data_particles[i]);
    free(phase_at_data_particles_perturb[i]);
  }
  free(Trans2D_at_veldata_particles);
  free(Trans2D_at_veldata_particles_perturb);
  free(TransTau_particles);
  free(TransTau_particles_perturb);
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);
  free(Trans_alpha_at_veldata_particles);
  free(Trans_alpha_at_veldata_particles_perturb);
  free(Trans_beta_at_veldata_particles);
  free(Trans_beta_at_veldata_particles_perturb);
  free(phase_at_data_particles);
  free(phase_at_data_particles_perturb);

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
  free(best_model_sarm);
  free(best_model_std_sarm);
  
  free(clouds_weight);
  free(clouds_alpha);
  free(clouds_beta);
  free(clouds_tau);
  free(clouds_vel);

  free(momentum_sarm_alpha);
  free(momentum_sarm_beta);

  free(workspace_phase);

  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  if(thistask == roottask)
  {
    printf("Ends reconstruct_sarm.\n");
  }

  return;
}

double prob_sarm(const void *model)
{
  double prob_sa = 0.0, prob_line=0.0, var2, dy;
  int i, j, k;
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
   * Trans2D is a pointer to the transfer function
   */
  if( (which_parameter_update < num_params_blr_tot) || force_update == 1 )
  {
    /* re-point */
    TransTau = TransTau_particles_perturb[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles_perturb[which_particle_update];
    Trans_alpha_at_veldata = Trans_alpha_at_veldata_particles_perturb[which_particle_update];
    Trans_beta_at_veldata = Trans_beta_at_veldata_particles_perturb[which_particle_update];

    transfun_sarm_cal_cloud((void *)pm, vel_sa_data_ext, Trans2D_at_veldata, Trans_alpha_at_veldata,
                            Trans_beta_at_veldata, n_vel_sarm_data_ext, 0);
  }
  else /* repoint */
  {
    TransTau = TransTau_particles[which_particle_update];
    Trans2D_at_veldata = Trans2D_at_veldata_particles[which_particle_update];
    Trans_alpha_at_veldata = Trans_alpha_at_veldata_particles[which_particle_update];
    Trans_beta_at_veldata = Trans_beta_at_veldata_particles[which_particle_update];
  }

  Fline2d_sarm_at_data = Fline_at_data_particles[which_particle_update];
  phase_sarm_at_data = phase_at_data_particles[which_particle_update];
  calculate_sarm_with_sample((void *)pm);

  /* loop over line */
  for(i=0; i<n_epoch_sarm_data; i++)
  {
    for(j=0; j<n_vel_sarm_data; j++)
    {
      //note mask with error < 0.0
      if(Flerrs2d_sarm_data[i*n_vel_sarm_data+j] > 0.0)
      {
        dy = Fline2d_sarm_data[i*n_vel_sarm_data + j] - Fline2d_sarm_at_data[i * n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)];
        var2 = Flerrs2d_sarm_data[i*n_vel_sarm_data+j]*Flerrs2d_sarm_data[i*n_vel_sarm_data+j];
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }

  /* loop over phase */
  for(i=0; i<n_epoch_sarm_data; i++)
  {
    for(k=0; k<n_base_sarm_data; k++)
    {
      for(j=0; j<n_vel_sarm_data; j++)
      {
        if(pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j] > 0.0)
        {
          dy = phase_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j] 
              -phase_sarm_at_data[i*n_vel_sarm_data_ext*n_base_sarm_data + k*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)];
          var2 = pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j]
                *pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j];
          prob_sa += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
  }
  return prob_sa + prob_line;
}

double prob_initial_sarm(const void *model)
{
  double prob_sa = 0.0, prob_line=0.0, var2, dy;
  int i, j, k;
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
  Trans_alpha_at_veldata = Trans_alpha_at_veldata_particles[which_particle_update];
  Trans_beta_at_veldata = Trans_beta_at_veldata_particles[which_particle_update];

  transfun_sarm_cal_cloud((void *)pm, vel_sa_data_ext, Trans2D_at_veldata, Trans_alpha_at_veldata,
                            Trans_beta_at_veldata, n_vel_sarm_data_ext, 0);
  
  Fline2d_sarm_at_data = Fline_at_data_particles[which_particle_update];
  phase_sarm_at_data = phase_at_data_particles[which_particle_update];
  calculate_sarm_with_sample((void *)pm);

  /* loop over line */
  for(i=0; i<n_epoch_sarm_data; i++)
  {
    for(j=0; j<n_vel_sarm_data; j++)
    {
      //note mask with error < 0.0
      if(Flerrs2d_sarm_data[i*n_vel_sarm_data+j] > 0.0)
      {
        /* note the vel extension */
        dy = Fline2d_sarm_data[i*n_vel_sarm_data + j] - Fline2d_sarm_at_data[i * n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)];
        var2 = Flerrs2d_sarm_data[i*n_vel_sarm_data+j]*Flerrs2d_sarm_data[i*n_vel_sarm_data+j];
        prob_line += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
      }
    }
  }

  /* loop over phase */
  for(i=0; i<n_epoch_sarm_data; i++)
  {
    for(k=0; k<n_base_sarm_data; k++)
    {
      for(j=0; j<n_vel_sarm_data; j++)
      {
        if(pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j] > 0.0)
        {
          /* note the vel extension */
          dy = phase_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j] 
              -phase_sarm_at_data[i*n_vel_sarm_data_ext*n_base_sarm_data + k*n_vel_sarm_data_ext + (j+n_vel_sa_data_incr)];
          var2 = pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j]
                *pherrs_sarm_data[i*n_vel_sarm_data*n_base_sarm_data + k*n_vel_sarm_data + j];
          prob_sa += (-0.5 * (dy*dy)/var2) - 0.5*log(var2 * 2.0*PI);
        }
      }
    }
  }

  return prob_sa + prob_line;
}
#endif