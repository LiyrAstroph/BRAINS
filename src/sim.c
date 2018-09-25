/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */


/*!
 * \file sim.c
 * \brief generate mocked 2d data.
 */
#include <stdio.h>
#include <stdlib.h> 
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

void *model;

void sim()
{
  if(thistask != roottask)
    return;

  FILE *fp;
  char fname[200];
  int i, j;

  sim_init();
  
  double *pm = (double *)model, error, fcon;
  
  switch(parset.flag_blrmodel)
  {
    case 1:
      pm[0] = log(1.0);
      pm[1] = 0.0;
      pm[2] = log(4.0);
      pm[3] = 0.9;
      pm[4] = 0.2;
      pm[5] = cos(20.0/180.0*PI);
      pm[6] = 40.0;
      pm[7] = 0.0;
      pm[8] = log(3.0);
      pm[9] = 0.1;
      pm[10] = 0.5;
      pm[11] = log(1.0);
      break;

    case 2:
      pm[0] = log(1.0);
      pm[1] = 0.0;
      pm[2] = log(4.0);
      pm[3] = 0.9;
      pm[4] = 0.2;
      pm[5] = cos(20.0/180.0*PI);
      pm[6] = 40.0;
      pm[7] = 0.0;
      pm[8] = log(3.0);
      pm[9] = log(0.01);
      pm[10] = log(0.1);
      pm[11] = log(1.0);
      break;

    case 3:
      pm[0] = log(1.0);
      pm[1] = 0.0;
      pm[2] = -1.0;
      pm[3] = log(3.0);
      pm[4] = log(5.0);
      pm[5] = cos(20.0/180.0*PI);
      pm[6] = 40.0;
      pm[7] = 0.0;
      pm[8] = log(3.0);
      pm[9] = 0.5;
      pm[10] = 0.5;
      pm[11] = log(1.0);

      break;

    case 4:
      pm[0] = log(1.0);
      pm[1] = 0.0;
      pm[2] = -1.0;
      pm[3] = log(3.0);
      pm[4] = log(5.0);
      pm[5] = cos(20.0/180.0*PI);
      pm[6] = 40.0;
      pm[7] = 0.0;
      pm[8] = log(3.0);
      pm[9] = 0.5;
      pm[10] = 0.5;
      pm[11] = log(1.0);
      break;

    case 5:
      pm[0] = log(1.0);  //A
      pm[1] = 0.0;       //Ag
      pm[2] = log(4.0);  //mu
      pm[3] = 0.5;       //Fin
      pm[4] = log(2.0);  //Fout
      pm[5] = 1.5;       //alpha
      pm[6] = cos(20.0/180.0*PI);      //inc
      pm[7] = 40.0;      //opn
      pm[8] = 0.5;       //k
      pm[9] = 2.0;       //gam
      pm[10] = 0.5;      //xi
      pm[11] = log(2.0); //mbh
      pm[12] = 0.5;      //fellip
      pm[13] = 0.5;      //fflow
      pm[14] = log(0.01);
      pm[15] = log(0.1);
      pm[16] = log(0.01);
      pm[17] = log(0.1);
      pm[18] = 0.0;      //theta_rot
      pm[19] = log(0.001); //sig_turb
      pm[20] = 0.0;      //spectral broadening
      pm[21] = log(1.0); //systematic error
      break;

    case 6:
      pm[0] = log(1.0);   // A
      pm[1] = 0.0;        // Ag
      pm[2] = log(4.0);   // mu
      pm[3] = 1.0;        // beta
      pm[4] = 0.25;        // F
      pm[5] = cos(20.0/180.0*PI); // inc
      pm[6] = 40.0;       // opn
      pm[7] = -0.4;        // kappa
      pm[8] = 5.0;        // gamma
      pm[9] = 0.5;        // obscuration
      pm[10] = log(2.0);  //mbh
      pm[11] = 0.5;       //fellip
      pm[12] = 0.4;       //fflow
      pm[13] = log(0.01); //
      pm[14] = log(0.1);  //
      pm[15] = log(0.01); //
      pm[16] = log(0.1);  //
      pm[17] = 0.0;       // theta_rot
      pm[18] = log(0.001);  // sig_turb
      pm[19] = 0.0;       // parameter for spectral broadening 
      pm[20] = log(1.0);  // systematic error
      break;

    case 7:
      pm[0] = log(1.0);  //A
      pm[1] = 0.0;       //Ag
      pm[2] = log(4.0);  //mu
      pm[3] = 0.8;       //beta
      pm[4] = 0.1;       //F
      pm[5] = cos(20.0/180.0*PI);//inc
      pm[6] = 40.0;      //opn
      pm[7] = 0.0;       //kappa
      pm[8] = 1.0;       //gamma
      pm[9] = 0.0;       //xi

      pm[10] = 0.5;      //fsh
      pm[11] = log(8.0); //mu_un
      pm[12] = 1.2;      //beta_un
      pm[13] = 0.1;      //F_un
      pm[14] = 20.0;     //opn_un

      pm[15] = log(2.0); //mbh
      pm[16] = 0.5;      //fellip
      pm[17] = 0.4;      //fflow
      pm[18] = log(0.01); 
      pm[19] = log(0.1);
      pm[20] = log(0.01);
      pm[21] = log(0.1);
      pm[22] = 0.0;      //theta_rot
      pm[23] = 0.5;      //fellip_un
      pm[24] = 0.4;      //fflow_un
      pm[25] = log(0.001); // sig_turb
      pm[26] = 0.0;      // parameter for spectral broadening 
      pm[27] = log(1.0); //systematic error
      break;
  }

  smooth_init(parset.n_vel_recon, TransV);
  
  if(parset.flag_dim == -1)
  {
    //note that here use sigma_hat = sigma/sqrt(tau).
    reconstruct_con_from_varmodel(exp(-3.714542), exp(3.104897), 1.0, 0.0); 
  }
  else
  {
    con_scale = 1.0;
    line_scale = 1.0;
    line_error_mean = con_error_mean = 0.01;
    create_con_from_random(0.03, 45.0, 1.0, 0.0); 
  }
  gsl_interp_init(gsl_linear, Tcon, Fcon, parset.n_con_recon);
  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_con_full.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i]/con_scale, Fcerrs[i]/con_scale);
  }
  fclose(fp);

  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_con.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  if(parset.flag_dim != -2)
  {
    for(i=0; i<n_con_data; i++)
    {
      //fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i]/con_scale, Fcerrs[i]/con_scale);
      fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, Tcon_data[i], gsl_acc);
      fprintf(fp, "%f %f %f\n", Tcon_data[i], (fcon+gsl_ran_ugaussian(gsl_r)*con_error_mean)/con_scale, con_error_mean/con_scale);
    }
  }
  else
  {
    for(i=0; i<parset.n_con_recon; i++)
    {
      fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i]/con_scale, Fcerrs[i]/con_scale);
    }
  }
  fclose(fp);
  
  
  transfun_1d_cal(model, 0);
  calculate_line_from_blrmodel(model, Tline, Fline, parset.n_line_recon);

  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_hb.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  if(parset.flag_dim == -1)
  {
    error = line_error_mean * sqrt(n_line_data) * (Vline_data[1] - Vline_data[0]);
  }
  else
  {
    error = line_error_mean;
  }
  for(i=0; i<parset.n_line_recon; i++)
  {
    fprintf(fp, "%f %f %f\n", Tline[i], Fline[i]/line_scale + gsl_ran_ugaussian(gsl_r)*error/line_scale, 
      error/line_scale);
  }

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
    fprintf(fp, "%f %f\n", TransTau[i], Trans1D[i]);
  }
  fclose(fp);

  transfun_2d_cal(model, TransV, Trans2D, parset.n_vel_recon, 1);
  calculate_line2d_from_blrmodel(model, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_hb2d.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  fprintf(fp, "# %d %d\n", parset.n_line_recon, parset.n_vel_recon);
  for(i=0; i<parset.n_line_recon; i++)
  {
    for(j=0; j<parset.n_vel_recon; j++)
    {
      fprintf(fp, "%f %f %f %f\n", TransV[j]*VelUnit, Tline[i],  
        (Fline2d[i*parset.n_vel_recon + j] + gsl_ran_ugaussian(gsl_r)*line_error_mean*0.3)/line_scale, line_error_mean/line_scale);
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
      fprintf(fp, "%f %f %f\n", TransV[j]*VelUnit, TransTau[i], Trans2D[i*parset.n_vel_recon + j]);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);

  smooth_end();
  sim_end();
}

void sim_init()
{
  int i;
  double dT, Tspan;

  switch(parset.flag_blrmodel)
  {
    case 1:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model1;
      break;
    case 2:
      num_params_blr_model = sizeof(BLRmodel2)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model2;
      break;
    case 3:
      num_params_blr_model = sizeof(BLRmodel3)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model3;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model3;
      break;
    case 4:
      num_params_blr_model = sizeof(BLRmodel4)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model3;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model4;
      break;
    case 5:
      num_params_blr_model = sizeof(BLRmodel5)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model5;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model5;
      break;

    case 6:
      num_params_blr_model = sizeof(BLRmodel6)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model6;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model6;
      break;

    case 7:
      num_params_blr_model = sizeof(BLRmodel7)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model7;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model7;
      break;

    default:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      transfun_1d_cloud_sample = transfun_1d_cloud_sample_model1;
      transfun_2d_cloud_sample = transfun_2d_cloud_sample_model1;
      break;
  }

  if(parset.InstRes <= 0.0)
  {
    num_params_res = 1;
    parset.InstRes = 220.0/VelUnit;
    parset.InstRes_err = 0.0;
  }

  if(parset.flag_narrowline > 1)
  {
    printf("# set flag_narrowline to 1.\n");
    parset.flag_narrowline = 1;
  }

  if(parset.flag_narrowline == 0)
  {
    parset.width_narrowline = 0.0;
  }

  parset.flag_linecenter = 0;

  parset.num_particles = 1;
  which_particle_update = 0;
  force_update = 1;
  which_parameter_update = -1;
  
  num_params_blr = num_params_blr_model + num_params_nlr + num_params_res + num_params_linecenter;
  num_params_var = 4 + parset.flag_trend;
  num_params = num_params_blr + num_params_var + parset.n_con_recon;

  model = malloc(num_params * sizeof(double));

  Fcon = malloc(parset.n_con_recon * sizeof(double));

  if(parset.flag_dim == -1)
  {
    Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  
  /* set time array for continuum */
    if(parset.time_back > 0.0)
      Tcon_min = Tcon_data[0] - parset.time_back;
    else
      Tcon_min = Tcon_data[0] - 0.5*Tspan;
    
    Tcon_max = Tcon_data[n_con_data-1] + fmax(0.05*Tspan, 10.0);
    dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
    for(i=0; i<parset.n_con_recon; i++)
    {
      Tcon[i] = Tcon_min + i*dT;
    }
  }
  else
  {
    Tspan = 120.0;

    rcloud_min_set = 0.0;
    rcloud_max_set = Tspan/2.0;
    if(parset.rcloud_max > 0.0)
      rcloud_max_set = fmin(rcloud_max_set, parset.rcloud_max);
    printf("rcloud_min_max_set: %f %f\n", rcloud_min_set, rcloud_max_set);

    Tcon_min = 0.0;
    Tcon_max = Tspan;

    dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
    for(i=0; i<parset.n_con_recon; i++)
    {
      Tcon[i] = Tcon_min + i*dT;
    } 
  }

  if(parset.flag_dim == -1)
  {
    parset.n_line_recon = n_line_data;
    parset.n_vel_recon = n_vel_data;
  }

  TransTau = malloc(parset.n_tau * sizeof(double));
  TransV = malloc(parset.n_vel_recon * sizeof(double));
  Trans1D = malloc(parset.n_tau * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  if(parset.flag_dim == -1)
  {
    memcpy(Tline, Tline_data, n_line_data * sizeof(double));
    memcpy(TransV, Vline_data, n_vel_data * sizeof(double));
  }
  else
  {
    Tline_min = Tcon_min + 30.0;
    Tline_max = Tcon_max;

    dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

    for(i=0; i<parset.n_line_recon; i++)
    {
      Tline[i] = Tline_min + i*dT;
    }
    Tline[parset.n_line_recon - 1] = Tline_max-1.0e-10;

    double vel_max_set, vel_min_set;
    vel_max_set = 3000.0/VelUnit;
    vel_min_set = - vel_max_set;
    double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

    for(i=0; i<parset.n_vel_recon; i++)
    {
      TransV[i] = vel_min_set + dVel*i;
    }
  }
  
  clouds_particles = malloc(parset.num_particles * sizeof(double *));
  clouds_particles_perturb = malloc(parset.num_particles * sizeof(double *));
  for(i=0; i<parset.num_particles; i++)
  {
    clouds_particles[i] = malloc(parset.n_cloud_per_task * sizeof(double));
    clouds_particles_perturb[i] = malloc(parset.n_cloud_per_task * sizeof(double));
  }

  con_q = malloc(nq * sizeof(double));

  tmp_tau = malloc(parset.n_cloud_per_task * sizeof(double));
  tmp_weight = malloc(parset.n_cloud_per_task * sizeof(double));
  tmp_vel = malloc(parset.n_cloud_per_task * parset.n_vel_per_cloud * sizeof(double));


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

  return;
}

void sim_end()
{
  int i;
  free(model);
  free(Fcon);

  free(TransTau);
  free(TransV);
  free(Tline);
  free(Fline);
  free(Fline2d);
  free(Trans2D);
  free(Trans1D);
  for(i=0; i<parset.num_particles; i++)
  {
    free(clouds_particles[i]);
    free(clouds_particles_perturb[i]);
  }
  free(clouds_particles);
  free(clouds_particles_perturb);

  free(con_q);

  free(tmp_tau);
  free(tmp_weight);
  free(tmp_vel);

  if(parset.flag_save_clouds && thistask == roottask)
  {
    fclose(fcloud_out);
  }
}
