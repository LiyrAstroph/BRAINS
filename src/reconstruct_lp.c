/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file reconstruct_lp.c
 *  \brief reconstruct line profile and recover BLR model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_interp.h>


#include "brains.h"

void *best_model_lp;      /*!< best model */
void *best_model_std_lp;  /*!< standard deviation of the best model */

/*!
 * postprocessing.
 */
void postprocesslp()
{
  char posterior_sample_file[BRAINS_MAX_STR_LENGTH];
  char posterior_sample_file_info[BRAINS_MAX_STR_LENGTH];
  void *posterior_sample, *posterior_sample_info, *post_model;
  double *pm, *pmstd;
  int size_of_modeltype = num_params * sizeof(double);

  best_model_lp = malloc(size_of_modeltype);
  best_model_std_lp = malloc(size_of_modeltype);

  if(thistask == roottask)
  {
    smooth_init(n_vel_data_ext, Vline_data_ext);
    char fname[200], buf[200];
    int i, j, num_ps;
    FILE *fp, *fpinfo, *flp;
    double info, info_max=-DBL_MAX;
    int idx_max;
    double *pinfo;

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

    strcpy(posterior_sample_file_info, "data/posterior_sample_info_lp.txt");
    fpinfo = fopen(posterior_sample_file_info, "r");
    if(fpinfo == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s.\n", posterior_sample_file_info);
      exit(0);
    }
    fgets(buf, 200, fpinfo);

    //file for line reconstruction
    sprintf(fname, "%s/%s", parset.file_dir, "data/lineprofile_rec.txt");
    flp = fopen(fname, "w");
    if(flp == NULL)
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
    posterior_sample_info = malloc(num_ps * sizeof(double));
    pinfo = (double *)posterior_sample_info;
    
    Fline_at_data = Fline_at_data_particles[0];
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

      if(fscanf(fpinfo, "%lf\n", (double *)posterior_sample_info+i) < 1)
      {
        fprintf(stderr, "# Error: Cannot read file %s.\n", posterior_sample_file_info);
        exit(0);
      }

      if(info_max < pinfo[i])
      {
        info_max = pinfo[i];
        idx_max = i;
      }

      //store model
      memcpy(posterior_sample+i*size_of_modeltype, post_model, size_of_modeltype);

      gen_cloud_sample(post_model, 0, 0);
      cal_line_profile_with_sample(post_model, Vline_data_ext, Fline_at_data, n_vel_data_ext);

      // output line profiles 
      for(j = 0; j<n_vel_data; j++)
      {
        fprintf(flp, "%e %e\n", Wline_data[j], Fline_at_data[j + n_vel_data_incr]/profile_scale);
      }
      fprintf(flp, "\n");
    }

    pm = (double *)best_model_lp;
    pmstd = (double *)best_model_std_lp;
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
    
    printf("Best params from mean and std of posterior sample:\n");
    for(j = 0; j<num_params_blr; j++)
      printf("%3d %10.6f +- %10.6f\n", j, *((double *)best_model_lp + j), 
                                             *((double *)best_model_std_lp+j) ); 
    
    printf("Best params from the highest likelihood:\n");
    memcpy(best_model_lp, (void *)((double *)posterior_sample + idx_max*num_params), size_of_modeltype);
    for(j = 0; j<num_params_blr; j++)
      printf("%3d %10.6f\n", j, *((double *)best_model_lp + j)); 

    fclose(fp);
    fclose(fpinfo);
    fclose(flp);
    free(post_model);
    free(posterior_sample);
    smooth_end(); 
  }
}

/*!
 * this function run dnest sampleing, reconstruct line profile using the best estimates for parameters.
 */
void reconstruct_lp()
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
  strcat(argv[argc++], "/data/restartlp_dnest.txt");

  if(parset.flag_restart == 1)
  {
    strcpy(argv[argc++], "-r");
    strcpy(argv[argc], parset.file_dir);
    strcat(argv[argc], "/");
    strcat(argv[argc++], "data/restartlp_dnest.txt");
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
  
  // sample tag 
  strcpy(argv[argc++], "-g");
  strcpy(argv[argc++], "_lp");
  
  reconstruct_lp_init();

  smooth_init(n_vel_data_ext, Vline_data_ext);
  dnest_lp(argc, argv);
  smooth_end();

  if(parset.flag_exam_prior != 1 && parset.flag_para_name != 1)
  {
    postprocesslp();

    // calculate light curves using the best model
    if(thistask == roottask)
    {
      /* create a wavelength grid */
      double *fline, *vel, *wave, dV, vmax;
      vel = malloc(parset.n_vel_recon*sizeof(double));
      wave = malloc(parset.n_vel_recon*sizeof(double));
      fline = malloc(parset.n_vel_recon*sizeof(double));
      
      vmax = fmax(fabs(Vline_data_ext[n_vel_data_ext-1]), fabs(Vline_data_ext[0]));
      dV = (vmax*2)/(parset.n_vel_recon+1);
      for(i=0; i<parset.n_vel_recon; i++)
      {
        vel[i] = -vmax + i*dV;
        wave[i] = (1.0 + vel[i]/C_Unit) * parset.linecenter * (1.0+parset.redshift);
      }

      smooth_init(parset.n_vel_recon, vel);
      FILE *flp;
      char fname[200];
      
      gen_cloud_sample(best_model_lp, 0, parset.flag_save_clouds);
      cal_line_profile_with_sample(best_model_lp, vel, fline, parset.n_vel_recon);
      
      //file for line reconstruction
      sprintf(fname, "%s/%s", parset.file_dir, parset.plineprofile_out_file);
      flp = fopen(fname, "w");
      if(flp == NULL)
      {
        fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
        exit(0);
      }

      // output line profiles 
      for(i = 0; i<parset.n_vel_recon; i++)
      {
        fprintf(flp, "%e %e\n", wave[i], fline[i]/profile_scale);
      }
      
      fclose(flp);
      smooth_end();
      free(vel);
      free(wave);
      free(Fline);
    }
  }

  reconstruct_lp_end();

  //clear up argv
  for(i=0; i<11; i++)
  {
    free(argv[i]);
  }
  free(argv);

  return;
}

void reconstruct_lp_init()
{
  int i;

  sprintf(dnest_options_file, "%s/%s", parset.file_dir, "param/OPTIONSLP");
  if(thistask == roottask)
  {
    get_num_particles(dnest_options_file);
  }
  MPI_Bcast(&parset.num_particles, 1, MPI_INT, roottask, MPI_COMM_WORLD);

  Fline_at_data_particles = malloc(parset.num_particles * sizeof(double *));
  Fline_at_data_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    Fline_at_data_particles[i] = malloc(n_vel_data_ext * sizeof(double));
    Fline_at_data_particles_perturb[i] = malloc(n_vel_data_ext * sizeof(double));
  }

  clouds_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));

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

void reconstruct_lp_end()
{
  int i;
  
  free(par_fix);
  free(par_fix_val);

  free(best_model_lp);
  free(best_model_std_lp);

  for(i=0; i<parset.num_particles; i++)
  {
    free(Fline_at_data_particles[i]);
    free(Fline_at_data_particles_perturb[i]);
  }
  free(Fline_at_data_particles);
  free(Fline_at_data_particles_perturb);
  
  for(i=0; i<num_params; i++)
  {
    free(par_range_model[i]);
    free(par_prior_gaussian[i]);
  }
  free(par_range_model);
  free(par_prior_gaussian);
  free(par_prior_model);
  
  free(clouds_weight);
  free(clouds_vel);

  if(parset.flag_save_clouds && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}

double prob_lp(const void *model)
{
  int i;
  double prob = 0.0, dy, var2;

  which_particle_update = dnest_get_which_particle_update();
  Fline_at_data = Fline_at_data_particles[which_particle_update];
  gen_cloud_sample(model, 0, 0);
  cal_line_profile_with_sample(model, Vline_data_ext, Fline_at_data, n_vel_data_ext);
  for(i=0; i<n_vel_data; i++)
  {
    if(Flerrs_data[i] > 0.0)
    {
      dy = Fline_data[i] - Fline_at_data[i+n_vel_data_incr];  /* note n_vel_data_ext */
      var2 = Flerrs_data[i]*Flerrs_data[i];
      prob += (-0.5 *  (dy*dy)/var2 - 0.5*log(var2 * 2.0*PI));
    }
  }
  return prob;
}

void cal_line_profile_with_sample(const void *pm, double *vel, double *fv, int nvel)
{
  int i, j, idV;
  double scale, norm;
  double *model=(double *)pm;
  scale = exp(model[num_params - 1]); /* the last one is scale factor */

  double dV, V, V_offset, flux_sum;
  dV = vel[1] - vel[0];
  
  for(i=0; i<nvel; i++)
  {
    fv[i] = 0.0;
  }
  
  norm = 0.0;
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    norm += clouds_weight[i] * parset.n_vel_per_cloud;
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i];
            
      V_offset = V + bin_offset * dV; /* bin type: center or left edge */
      if(V_offset < vel[0] || V_offset >= vel[nvel-1] + dV )
        continue;
      idV = (V_offset - vel[0])/dV; 
      fv[idV] += clouds_weight[i];
    }
  }
  norm *= dV;

  for(i=0; i<nvel; i++)
  {
    fv[i] = fv[i] / norm * scale;
  }

  /* add intrinsic narrow line */
  if(parset.flag_narrowline != 0)
  {
    double flux, width, shift, fnarrow;
    if(parset.flag_narrowline == 1)  /* fixed narrow line */
    {
      flux = parset.flux_narrowline;
      width = parset.width_narrowline;
      shift = parset.shift_narrowline;
    }
    else if(parset.flag_narrowline == 2) /* narrow line with Gaussian priors */
    {
      flux =  parset.flux_narrowline  + model[num_params_blr_model] * parset.flux_narrowline_err;
      width = parset.width_narrowline + model[num_params_blr_model+1] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + model[num_params_blr_model+2] * parset.shift_narrowline_err;
    }
    else  /* narrow line with logrithmic prior of flux */
    {
      flux =  exp(model[num_params_blr_model]);
      width = parset.width_narrowline + model[num_params_blr_model+1] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + model[num_params_blr_model+2] * parset.shift_narrowline_err;
    }

    width = fmax(1.0e-10, width); /* make sure thant width is not zero */

    for(i=0; i<nvel; i++)
    {
      fnarrow = flux/sqrt(2.0*PI)/width * exp( -0.5 * pow( (vel[i] - shift)/(width), 2.0) );
      fv[i] += fnarrow;
    } 
  }
  
  line_gaussian_smooth_FFT(vel, fv, nvel, pm);

  return;
}