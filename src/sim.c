/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
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
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>
#include <mpi.h>

#include "brains.h"

void *model;

void sim()
{
  if(thistask != roottask)
    return;

  FILE *fp;
  char fname[200];
  int i, j, k, incr, flag_save_clouds=0;

  sim_init();
  
  double *pm = (double *)model, error, fcon;
  double sigma, taud;
  
  smooth_init(parset.n_vel_recon, TransV);

/*==================================continuum============================================*/  
  if(parset.flag_dim == -1)
  {
    /* note that here use sigma_hat = sigma/sqrt(tau) */
    sigma = var_param[1];
    taud = var_param[2];
    printf("Sim with ln(sigma) = %f and  ln(taud) = %f.\n", sigma, taud);
    reconstruct_con_from_varmodel(exp(var_param[1]), exp(var_param[2]), 1.0, 0.0); 
  }
  else
  {
    con_scale = 1.0;
    line_scale = 1.0;
    line_error_mean = con_error_mean = 0.01;
    /* arguments: sigma_hat, tau, alapha, and syserr 
     * note the light curve is shifted and scaled, so that sigma_hat maybe changed. */
    sigma = 0.03;
    taud = (Tcon[parset.n_con_recon-1] - Tcon[0])/10.0;
    printf("Sim with sigma = %f and  taud = %f.\n", sigma, taud);
    create_con_from_random(sigma, taud, 1.0, 0.0); 
  }
  calculate_con_rm(model);
  gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);
  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_con_full.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  for(i=0; i<parset.n_con_recon; i++)
  {
    fprintf(fp, "%e %e %e\n", Tcon[i]*(1.0+parset.redshift), Fcon[i]/con_scale, Fcerrs[i]/con_scale);
  }
  fclose(fp);

  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_con.txt");
  printf("%-13s %s\n", "Continuum:", fname);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }
  
  if(parset.flag_dim != -2)
  {
    gsl_interp_init(gsl_linear, Tcon, Fcon_rm, parset.n_con_recon);
    for(i=0; i<n_con_data; i++)
    {
      //fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i]/con_scale, Fcerrs[i]/con_scale);
      fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, Tcon_data[i], gsl_acc);
      fprintf(fp, "%e %e %e\n", Tcon_data[i]*(1.0+parset.redshift), 
            (fcon+gsl_ran_ugaussian(gsl_r)*con_error_mean)/con_scale, con_error_mean/con_scale);
    }
  }
  else
  {
    incr = fmax(0.5/(Tcon[1]-Tcon[0]), 1.0); //cadence to be 0.5day or increasement to be 1
    for(i=0; i<parset.n_con_recon; i+=incr)
    {
      if(Tcon[i] >= 0.0)
      {
        fprintf(fp, "%e %e %e\n", Tcon[i]*(1.0+parset.redshift), 
                (Fcon[i]+gsl_ran_ugaussian(gsl_r)*Fcerrs[i])/con_scale, Fcerrs[i]/con_scale);
      }
    }
  }
  fclose(fp);

/*==================================1d line============================================*/  
  transfun_1d_cal(model, flag_save_clouds);
  calculate_line_from_blrmodel(model, Tline, Fline, parset.n_line_recon);
  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_line.txt");
  printf("%-13s %s\n", "Line:", fname);
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
    fprintf(fp, "%e %e %e\n", Tline[i]*(1.0+parset.redshift), Fline[i]/line_scale + gsl_ran_ugaussian(gsl_r)*error/line_scale, 
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
    fprintf(fp, "%e %e\n", TransTau[i], Trans1D[i]);
  }
  fclose(fp);

/*==================================2d line============================================*/  
  flag_save_clouds = 0;
#ifndef SpecAstro
  if(parset.flag_blrmodel == -1) /* in this case, line profile sim does not produce clouds */
  {
    flag_save_clouds = 1;
  }
#endif
  transfun_2d_cal(model, TransV, Trans2D, parset.n_vel_recon, flag_save_clouds);
  calculate_line2d_from_blrmodel(model, Tline, TransV, 
          Trans2D, Fline2d, parset.n_line_recon, parset.n_vel_recon);

  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_line2d.txt");
  printf("%-13s %s\n", "Line2d:", fname);
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
      fprintf(fp, "%e %e %e\n", TransW[j],  
        (Fline2d[i*parset.n_vel_recon + j] + gsl_ran_ugaussian(gsl_r)*line_error_mean*0.3)/line_scale, line_error_mean/line_scale);
    }

    fprintf(fp, "\n");
  }
  fclose(fp);
  
  sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_broadening.txt");
  printf("%-13s %s\n", "Line Broaden:", fname);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s\n", fname);
    exit(-1);
  }

  for(i=0; i<parset.n_line_recon; i++)
  {
    fprintf(fp, "%f %f\n", parset.InstRes * VelUnit, parset.InstRes_err * VelUnit);
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

/*==================================line profile============================================*/
  if(parset.flag_blrmodel != -1)  
  {
    flag_save_clouds = 0;
#ifndef SpecAstro
    flag_save_clouds = 1;
#endif
    gen_cloud_sample(model, 0, flag_save_clouds);
    cal_line_profile_with_sample(model, TransV, Fline2d, parset.n_vel_recon);
    sprintf(fname, "%s/%s", parset.file_dir, "/data/sim_lineprofile.txt");
    printf("%-13s %s\n", "Line Profile:", fname);
    fp = fopen(fname, "w");
    if(fp == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
    
    fprintf(fp, "# %d\n", parset.n_vel_recon);
    for(j=0; j<parset.n_vel_recon; j++)
    {
      fprintf(fp, "%e %e %e\n", TransW[j],  
        (Fline2d[j] + gsl_ran_ugaussian(gsl_r)*line_error_mean*0.3)/line_scale, line_error_mean/line_scale);
    }
    fclose(fp);
  }

#ifdef SpecAstro
  double *sa_pm;
  sa_pm = (double *)pm + num_params_blr;
/*==================================SA============================================*/  
  sa_smooth_init(parset.n_sa_vel_recon, vel_sa, parset.sa_InstRes);

  gen_sa_cloud_sample((void *)sa_pm, 3, flag_save_clouds);
  calculate_sa_sim_with_sample(pm, vel_sa, parset.n_sa_vel_recon, base_sa, parset.n_sa_base_recon, 
                                   phase_sa, Fline_sa);
  
  sprintf(fname, "%s/%s", parset.file_dir, "data/sim_sa.txt");
  printf("%-13s %s\n", "SA:", fname);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  // output sa line
  fprintf(fp, "# %d %d %d\n", 1, parset.n_sa_vel_recon, parset.n_sa_base_recon);
  for(j=0; j<parset.n_sa_vel_recon; j++)
  {
    fprintf(fp, "%e %e %e\n", wave_sa[j], Fline_sa[j] + gsl_ran_ugaussian(gsl_r)*sa_line_error_mean, sa_line_error_mean);
  }
  fprintf(fp, "\n");
  for(i=0; i<parset.n_sa_base_recon; i++)
  {
    fprintf(fp, "# %f %f\n", base_sa[i*2], base_sa[i*2+1]);
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e %e\n", wave_sa[j], 
       phase_sa[i*parset.n_sa_vel_recon + j]/(ScaleFactor[j]) + gsl_ran_ugaussian(gsl_r)*sa_phase_error_mean, 
       sa_phase_error_mean);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  sa_smooth_end();

/*==================================SARM============================================*/  
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    Fcon_sarm[i] = gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, Tline_sarm[i], gsl_acc);
  }
  
  flag_save_clouds = 1;
  sarm_smooth_init(parset.n_sa_vel_recon, vel_sa, parset.sa_InstRes);
  transfun_sarm_cal_cloud((void *)sa_pm, vel_sa, Trans2D, Trans_sarm_alpha, Trans_sarm_beta, parset.n_sa_vel_recon, flag_save_clouds);
  calculate_sarm_sim_with_sample(pm, Tline_sarm, vel_sa, Trans2D, Trans_sarm_alpha, Trans_sarm_beta, 
                                     parset.n_sa_vel_recon, parset.n_sarm_line_recon, base_sarm, 
                                     parset.n_sarm_base_recon, phase_sarm, Fline_sarm, 
                                     momentum_sarm_alpha, momentum_sarm_beta,
                                     photocenter_sarm_alpha, photocenter_sarm_beta); 

  sprintf(fname, "%s/%s", parset.file_dir, "data/sim_sarm.txt");
  printf("%-13s %s\n", "SARM:", fname);
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  
  fprintf(fp, "# %d %d %d\n", parset.n_sarm_line_recon, parset.n_sa_vel_recon, parset.n_sarm_base_recon);
  /* output line */
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    fprintf(fp, "# %f  %e\n", Tline_sarm[i]*(1.0+parset.redshift), Fcon_sarm[i] + gsl_ran_ugaussian(gsl_r) * con_error_mean);
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e %e\n", wave_sa[j], Fline_sarm[i*parset.n_sa_vel_recon + j] + gsl_ran_ugaussian(gsl_r)*sarm_line_error_mean, 
                                            sarm_line_error_mean);
    }
    fprintf(fp,"\n");
  }
  /* output phase */
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    for(k=0; k<parset.n_sarm_base_recon; k++)
    {
      fprintf(fp, "# %f %f\n", base_sarm[i*parset.n_sarm_base_recon*2+k*2+0], base_sarm[i*parset.n_sarm_base_recon*2+k*2+1]);
      for(j=0; j<parset.n_sa_vel_recon; j++)
      {
        fprintf(fp, "%e %e %e\n", wave_sa[j], 
        phase_sarm[i*parset.n_sa_vel_recon*parset.n_sarm_base_recon + k*parset.n_sa_vel_recon + j]/(ScaleFactor[j]) 
        + gsl_ran_ugaussian(gsl_r)*sarm_phase_error_mean,
        sarm_phase_error_mean);
      }
      fprintf(fp,"\n");
    }
  }
  fclose(fp);

  /* output SARM transfer function */
  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_tran_alpha.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e %e\n", vel_sa[j]*VelUnit, TransTau[i], Trans_sarm_alpha[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_tran_beta.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e %e\n", vel_sa[j]*VelUnit, TransTau[i], Trans_sarm_beta[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  /* output SARM photocenter function */
  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_photocenter_alpha.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e\n", vel_sa[j]*VelUnit, photocenter_sarm_alpha[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_photocenter_beta.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e\n", vel_sa[j]*VelUnit, photocenter_sarm_beta[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  /* output SARM photocenter function */
  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_momentum_alpha.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e\n", vel_sa[j]*VelUnit, momentum_sarm_alpha[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sprintf(fname, "%s/%s", parset.file_dir, "data/sarm_momentum_beta.txt");
  fp = fopen(fname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "# Error: Cannot open file %s.\n", fname);
    exit(0);
  }
  for(i=0; i<parset.n_sarm_line_recon; i++)
  {
    for(j=0; j<parset.n_sa_vel_recon; j++)
    {
      fprintf(fp, "%e %e\n", vel_sa[j]*VelUnit, momentum_sarm_beta[i*parset.n_sa_vel_recon + j]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);

  sarm_smooth_end();
  
#endif  

  smooth_end();
  sim_end();
}

void sim_init()
{
  int i, j, idx;
  double dT, Tspan;
  double *pm, Rblr, mbh;

  switch(parset.flag_blrmodel)
  {
    case -1:
      num_params_blr_model = num_params_MyTransfun2d;
      transfun_1d_cal = transfun_1d_cal_mytransfun;
      transfun_2d_cal = transfun_2d_cal_mytransfun;
      BLRmodel_name = MyTFmodel_name;
      break;

    case 0:
      num_params_blr_model = num_params_MyBLRmodel2d;
      gen_cloud_sample = gen_cloud_sample_mymodel;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = MyBLRmodel_name;
      break;

    case 1:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;
    case 2:
      num_params_blr_model = sizeof(BLRmodel2)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model2;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel2_name;
      break;
    case 3:
      num_params_blr_model = sizeof(BLRmodel3)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model3;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel3_name;
      break;
    case 4:
      num_params_blr_model = sizeof(BLRmodel4)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model4;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel4_name;
      break;
    case 5:
      num_params_blr_model = sizeof(BLRmodel5)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model5;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel5_name;
      break;

    case 6:
      num_params_blr_model = sizeof(BLRmodel6)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model6;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel6_name;
      break;

    case 7:
      num_params_blr_model = sizeof(BLRmodel7)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model7;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel7_name;
      break;

    case 8:
      num_params_blr_model = sizeof(BLRmodel8)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model8;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel8_name;
      break;
    
    case 9:
      num_params_blr_model = sizeof(BLRmodel9)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model9;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel9_name;
      break;

    default:
      num_params_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_cloud_sample = gen_cloud_sample_model1;
      transfun_1d_cal = transfun_1d_cal_cloud;
      transfun_2d_cal = transfun_2d_cal_cloud;
      BLRmodel_name = BLRmodel1_name;
      break;
  }

#ifdef SpecAstro
  set_sa_blr_model();
  /* SA */
  num_params_sa_blr = num_params_sa_blr_model + num_params_sa_extpar;
  num_params_sa = num_params_sa_blr;

#endif
  
  /* use epoch-independent broadening */
  if(parset.flag_InstRes > 2)
  {
    num_params_res = 1;
    parset.InstRes = 220.0/VelUnit;
    parset.InstRes_err = 0.0;
  }
  else
  {
    parset.InstRes /= VelUnit;
    parset.InstRes_err = 0.0;
  }

  if(parset.flag_narrowline > 1)
  {
    printf("# set flag_narrowline to 1.\n");
    parset.flag_narrowline = 1;
  }
  else if(parset.flag_narrowline == 0)
  {
    parset.width_narrowline = 0.0;
  }

  parset.flag_linecenter = 0;
  num_params_linecenter = 0;

  parset.num_particles = 1;
  which_particle_update = 0;
  force_update = 1;
  which_parameter_update = -1;
  
  num_params_blr = num_params_blr_model + num_params_nlr 
                 + num_params_res + num_params_linecenter + 1; /* include line sys err */
  num_params_var = num_params_drw + num_params_trend + num_params_difftrend + num_params_resp;

  num_params_blr_tot = num_params_blr;
#ifdef SpecAstro
  num_params_blr_tot  += num_params_sa_blr;
#endif
  num_params = num_params_blr_tot + num_params_var + parset.n_con_recon;

  /* index of A and Ag */
  idx_resp = num_params_blr_tot + num_params_drw + num_params_trend;
  /* index of different trend */
  idx_difftrend = idx_resp + num_params_resp;

  model = malloc(num_params * sizeof(double));
  par_fix = (int *) malloc(num_params * sizeof(int));
  par_fix_val = (double *) malloc(num_params * sizeof(double));

  /* setup parameters */
  pm = (double *)model;
  if(parset.flag_blrmodel == -1)
  {
    set_par_value_mytransfun_sim(pm);
  }
  else 
  {
    if(parset.flag_para_value == 0)
    {
      set_par_value_sim(pm, parset.flag_blrmodel);
    }
    else 
    {
      read_param_value(pm, parset.para_value_file);
    }
    set_par_fix_blrmodel();
    for(i=0; i<num_params_blr_model; i++)
    {
      if(par_fix[i] == 1)
      {
        pm[i] = par_fix_val[i];
      }
    }
  }
  printf("=================================\n");
  printf("RM model parameter values:\n");
  print_par_value_sim(pm, num_params_blr_model);
  printf("=================================\n");
  
  /* spectral broadening, note this is a deviation from the input value */
  pm[num_params_blr_model + num_params_nlr ] = 0.0; 
  pm[num_params_blr -1] = 0.0; //line syserr=0 (for flag_dim!=3) or scale=1 (for flag_dim=3)

  pm[idx_resp + 0] = log(1.0); //A
  pm[idx_resp + 1] = 0.0;      //Ag

#ifdef SpecAstro
  double *sa_model = pm + num_params_blr;
  set_idx_par_mutual();
  if(parset.flag_para_value == 0)
  {
    set_par_value_sim(sa_model, parset.flag_sa_blrmodel);
  }
  else 
  {
    read_param_value(sa_model, parset.para_value_file);
  }
  set_par_fix_sa_blrmodel(); 
  for(i=0; i<num_params_sa_blr_model; i++)
  {
    if(par_fix[num_params_blr + i] == 1)
    {
      pm[num_params_blr + i] = par_fix_val[num_params_blr + i];
    }
  }

  /* set the same mbh and inc */
  sa_model[idx_sa_par_mutual[0]] = pm[idx_rm_par_mutual[0]]; //mbh
  sa_model[idx_sa_par_mutual[1]] = pm[idx_rm_par_mutual[1]]; //inc

  sa_model[num_params_sa_blr_model]     = log(550.0); //DA
  sa_model[num_params_sa_blr_model + 1] = 0.0;        //PA
  sa_model[num_params_sa_blr_model + 2] = 0.0;        //FA
  sa_model[num_params_sa_blr_model + 3] = 0.0;        //CO
  sa_model[num_params_sa_blr_model + 4] = 0.0;        //alphac
  sa_model[num_params_sa_blr_model + 5] = 0.0;        //betac
#endif

  Fcon = malloc(parset.n_con_recon * sizeof(double));
  Fcon_rm = malloc(parset.n_con_recon * sizeof(double));

  idx = get_idx_mbh_from_blrmodel(parset.flag_blrmodel);
  mbh = exp(pm[idx]);

  if(parset.flag_blrmodel != 8) /* model 8 is particular */
  {
    Rblr = exp(pm[0]);
  }
  else 
  {
    Rblr = exp(pm[9]);
  }

  if(parset.flag_dim == -1)
  {
    Tspan = Tcon_data[n_con_data -1] - Tcon_data[0];
  
  /* set time array for continuum */
    Tcon_min = Tcon_data[0] - time_back_set - 10.0;
    
    Tcon_max = fmax(Tcon_data[n_con_data-1], Tline_data[n_line_data-1]) + fmax(0.05*Tspan, 10.0);
    dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
    for(i=0; i<parset.n_con_recon; i++)
    {
      Tcon[i] = Tcon_min + i*dT;
    }
  }
  else
  {
    Tspan = Rblr*10.0;

    rcloud_min_set = 0.0;
    rcloud_max_set = Tspan/2.0;
    if(parset.rcloud_max > 0.0)
      rcloud_max_set = fmin(rcloud_max_set, parset.rcloud_max);
    printf("RM rcloud_min_max_set: %f %f\n", rcloud_min_set, rcloud_max_set);

    time_back_set = 2.0*rcloud_max_set;
    Tcon_min = 0.0 - time_back_set - 10.0;
    Tcon_max = Tspan + 10.0;

    dT = (Tcon_max - Tcon_min)/(parset.n_con_recon -1);
  
    for(i=0; i<parset.n_con_recon; i++)
    {
      Tcon[i] = Tcon_min + i*dT;
    } 
  }

  /* set Larr_rec */
  for(i=0;i<parset.n_con_recon;i++)
  {
    Larr_rec[i*nq + 0]=1.0;
    for(j=1; j<nq; j++)
      Larr_rec[i*nq + j] = pow(Tcon[i], j);
  }

  if(parset.flag_dim == -1)
  {
    parset.n_line_recon = n_line_data;
    parset.n_vel_recon = n_vel_data;
  }

  TransTau = malloc(parset.n_tau * sizeof(double));
  TransV = malloc(parset.n_vel_recon * sizeof(double));
  TransW = malloc(parset.n_vel_recon * sizeof(double));
  Trans1D = malloc(parset.n_tau * sizeof(double));
  Trans2D = malloc(parset.n_tau * parset.n_vel_recon * sizeof(double));
  Tline = malloc(parset.n_line_recon * sizeof(double));
  Fline = malloc(parset.n_line_recon * sizeof(double));
  Fline2d = malloc(parset.n_line_recon * parset.n_vel_recon * sizeof(double));

  if(parset.flag_dim == -1)
  {
    memcpy(Tline, Tline_data, n_line_data * sizeof(double));
    memcpy(TransV, Vline_data, n_vel_data * sizeof(double));
    memcpy(TransW, Wline_data, n_vel_data * sizeof(double));
  }
  else
  {
    Tline_min = 0.0;
    Tline_max = Tcon_max - 1.0;

    dT = (Tline_max - Tline_min)/(parset.n_line_recon - 1);

    for(i=0; i<parset.n_line_recon; i++)
    {
      Tline[i] = Tline_min + i*dT;
    }

    double vel_max_set, vel_min_set;
    vel_max_set = sqrt(pow(2.0*sqrt(mbh/Rblr), 2.0) + pow(2.0*parset.InstRes, 2.0));
    vel_min_set = - vel_max_set;
    double dVel = (vel_max_set- vel_min_set)/(parset.n_vel_recon -1.0);

    for(i=0; i<parset.n_vel_recon; i++)
    {
      TransV[i] = vel_min_set + dVel*i;
      TransW[i] = (1.0 + TransV[i]/C_Unit) * parset.linecenter * (1.0+parset.redshift);
    }
  }
  
  clouds_tau = malloc(parset.n_cloud_per_task * sizeof(double));
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

#ifdef SpecAstro
  double saRblr, norm;

  if(parset.flag_dim == -1)
  {
    parset.n_sa_vel_recon = n_vel_sa_data;
    parset.n_sa_base_recon = n_base_sa_data;

    parset.n_sarm_base_recon = n_base_sarm_data;
    parset.n_sarm_line_recon = n_epoch_sarm_data;
  }
  else
  {
    sa_flux_norm = 1.0;
    parset.n_sa_vel_recon = parset.n_vel_recon;  /* note SARM also uses tran2d, so no larger than "n_vel_recon" */
    if(parset.flag_gravity == 1)
      parset.n_sa_base_recon = n_base_sa_3c273;
    else
      parset.n_sa_base_recon = 20;

    sa_phase_error_mean = 0.01;
    sa_line_error_mean = 0.01;

    sarm_phase_error_mean = 0.01;
    sarm_line_error_mean = 0.01;

    parset.n_sarm_base_recon = 2;
    parset.n_sarm_line_recon = parset.n_line_recon;

    ScaleFactor = malloc(parset.n_sa_vel_recon * sizeof(double));
  }
    
  vel_sa = malloc(parset.n_sa_vel_recon * sizeof(double));
  wave_sa = malloc(parset.n_sa_vel_recon * sizeof(double));
  Fline_sa = malloc(parset.n_sa_vel_recon * sizeof(double));
  base_sa = malloc(parset.n_sa_base_recon * 2 * sizeof(double));
  phase_sa = malloc(parset.n_sa_vel_recon * parset.n_sa_base_recon * sizeof(double));
  
  clouds_alpha = malloc(parset.n_cloud_per_task * sizeof(double));
  clouds_beta = malloc(parset.n_cloud_per_task * sizeof(double));

  /* sarm */
  Tline_sarm = malloc(parset.n_sarm_line_recon * sizeof(double));
  Fcon_sarm = malloc(parset.n_sarm_line_recon * sizeof(double));
  Fline_sarm = malloc(parset.n_sarm_line_recon * parset.n_sa_vel_recon * sizeof(double));
  base_sarm = malloc(parset.n_sarm_base_recon * parset.n_sarm_line_recon * 2 * sizeof(double));
  momentum_sarm_alpha = malloc(parset.n_sa_vel_recon * parset.n_sarm_line_recon * sizeof(double));
  momentum_sarm_beta = malloc(parset.n_sa_vel_recon * parset.n_sarm_line_recon * sizeof(double));
  photocenter_sarm_alpha = malloc(parset.n_sa_vel_recon * parset.n_sarm_line_recon * sizeof(double));
  photocenter_sarm_beta = malloc(parset.n_sa_vel_recon * parset.n_sarm_line_recon * sizeof(double));
  phase_sarm = malloc(parset.n_sa_vel_recon * parset.n_sarm_base_recon * parset.n_sarm_line_recon * sizeof(double));
  Trans_sarm_alpha = malloc(parset.n_tau * parset.n_sa_vel_recon * sizeof(double));
  Trans_sarm_beta = malloc(parset.n_tau * parset.n_sa_vel_recon * sizeof(double));

  workspace_phase = malloc(parset.n_sa_vel_recon * 3 * sizeof(double));
  
  if(parset.flag_sa_blrmodel != 8)
  {
    saRblr = exp(pm[num_params_blr]);
  }
  else 
  {
    saRblr = exp(pm[num_params_blr + 9]);
  }
  rcloud_max_set = fmax(rcloud_max_set, saRblr * 5.0);
  printf("SA rcloud_min_max_set: %f %f\n", rcloud_min_set, rcloud_max_set);

  if(parset.flag_dim == -1)
  {
    /* SA */
    memcpy(vel_sa, vel_sa_data, parset.n_sa_vel_recon * sizeof(double));
    memcpy(wave_sa, wave_sa_data, parset.n_sa_vel_recon * sizeof(double));

    memcpy(base_sa, base_sa_data, parset.n_sa_base_recon * 2 * sizeof(double));
    
    /* SARM */
    memcpy(base_sarm, base_sarm_data, n_epoch_sarm_data*n_base_sarm_data*2 * sizeof(double));
    memcpy(Tline_sarm, Tline_sarm_data, n_epoch_sarm_data * sizeof(double));
  }
  else
  {
    double vel_max_set, vel_min_set, phi;
    vel_max_set = sqrt(pow(2.0*sqrt(mbh/saRblr), 2.0) + pow(2.0*parset.sa_InstRes, 2.0));
    vel_min_set = - vel_max_set;
    double dVel = (vel_max_set- vel_min_set)/(parset.n_sa_vel_recon -1.0);
    for(i=0; i<parset.n_sa_vel_recon; i++)
    {
      vel_sa[i] = vel_min_set + dVel*i;
      wave_sa[i] = (1.0 + vel_sa[i]/C_Unit) * parset.sa_linecenter * (1.0+parset.redshift);
    }
    
    if(parset.flag_gravity == 1)
      memcpy(base_sa, base_sa_3c273, n_base_sa_3c273*2*sizeof(double));
    else
    {
      for(i=0; i<parset.n_sa_base_recon; i++)
      {
        phi = -PI/2.0 + PI/parset.n_sa_base_recon * i;
        base_sa[i*2+0] = 100.0*cos(phi);
        base_sa[i*2+1] = 100.0*sin(phi);
      }
    }
    
    /* convert to direction */
    if(parset.flag_sa_datatype == 1)
    {
      for(i=0; i<parset.n_sa_base_recon; i++)
      {
        norm = pow(base_sa[i*2+0], 2) + pow(base_sa[i*2+1], 2);
        base_sa[i*2+0] /= sqrt(norm);
        base_sa[i*2+1] /= sqrt(norm);
      }
    }
    
    /* SARM */
    for(i=0; i<parset.n_sarm_line_recon; i++)
    {
      Tline_sarm[i] = Tline[i];
      for(j=0; j<parset.n_sarm_base_recon; j++)
      {
        base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 0] =  70.0 + 10*cos(PI/parset.n_sarm_base_recon * j);
        base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 1] = -70.0 + 10*sin(PI/parset.n_sarm_base_recon * j);
        
        /* convert to direction */
        if(parset.flag_sa_datatype == 1)
        {
          norm = pow(base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 0], 2)
                +pow(base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 1], 2);
          
          base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 0] /= sqrt(norm);
          base_sarm[i*parset.n_sarm_base_recon*2 + j*2 + 1] /= sqrt(norm);
        }
      }
    }

    /* setup scale factor */
    if(parset.flag_sa_datatype == 0)
    {
      sign = -1; /* phase */
      for(i=0; i<parset.n_sa_vel_recon; i++)
      {
        ScaleFactor[i] = PhaseFactor * wave_sa[i];
      }
    }
    else 
    {
      sign = 1; /* photocenter */
      for(i=0; i<parset.n_sa_vel_recon; i++)
      {
        ScaleFactor[i] = PhotoFactor;
      }
    }
  }
#endif

  return;
}

void sim_end()
{
  free(model);
  free(par_fix);
  free(par_fix_val);
  free(Fcon);
  free(Fcon_rm);

  free(TransTau);
  free(TransV);
  free(TransW);
  free(Tline);
  free(Fline);
  free(Fline2d);
  free(Trans2D);
  free(Trans1D);

  free(clouds_tau);
  free(clouds_weight);
  free(clouds_vel);

  if(parset.flag_save_clouds && thistask == roottask)
  {
    fclose(fcloud_out);
  }

#ifdef SpecAstro

  if(parset.flag_dim == -2)
    free(ScaleFactor);

  free(vel_sa);
  free(wave_sa);
  free(Fline_sa);
  free(base_sa);
  free(phase_sa);
  
  free(clouds_alpha);
  free(clouds_beta);
 
  free(Tline_sarm);
  free(Fcon_sarm);
  free(Fline_sarm);
  free(base_sarm);
  free(phase_sarm);
  free(momentum_sarm_alpha);
  free(momentum_sarm_beta);
  free(photocenter_sarm_alpha);
  free(photocenter_sarm_beta);
  free(Trans_sarm_alpha);
  free(Trans_sarm_beta);
  
  free(workspace_phase);
#endif
}

/* 
 * set parameter values for either RM or SA
 */
void set_par_value_sim(double *pm, int flag_model)
{
  int i;
  switch(flag_model)
  {
    case -1:
      set_par_value_mytransfun_sim(pm);
      break;

    case 0:
      set_par_value_mymodel_sim(pm);
      break;

    case 1:
      i=0;
      pm[i++] = log(4.0);
      pm[i++] = 0.9;
      pm[i++] = 0.2;
      pm[i++] = cos(20.0/180.0*PI);
      pm[i++] = 40.0;
      pm[i++] = 0.0;
      pm[i++] = log(3.0);
      pm[i++] = 0.1;
      pm[i++] = 0.5;
      break;

    case 2:
      i=0;
      pm[i++] = log(4.0);
      pm[i++] = 0.9;
      pm[i++] = 0.2;
      pm[i++] = cos(20.0/180.0*PI);
      pm[i++] = 40.0;
      pm[i++] = 0.0;
      pm[i++] = log(3.0);
      pm[i++] = log(0.01);
      pm[i++] = log(0.1);
      break;

    case 3:
      i=0;
      pm[i++] = log(3.0);
      pm[i++] = log(5.0);
      pm[i++] = -1.0;
      pm[i++] = cos(20.0/180.0*PI);
      pm[i++] = 40.0;
      pm[i++] = 0.0;
      pm[i++] = log(3.0);
      pm[i++] = 0.5;
      pm[i++] = 0.5;
      break;

    case 4:
      i=0;
      pm[i++] = log(3.0);
      pm[i++] = log(5.0);
      pm[i++] = -1.0;
      pm[i++] = cos(20.0/180.0*PI);
      pm[i++] = 40.0;
      pm[i++] = 0.0;
      pm[i++] = log(3.0);
      pm[i++] = 0.5;
      pm[i++] = 0.5;
      break;

    case 5:
      i=0;
      pm[i++] = log(4.0);  //mu
      pm[i++] = 0.5;       //Fin
      pm[i++] = log(2.0);  //Fout
      pm[i++] = 1.5;       //alpha
      pm[i++] = cos(20.0/180.0*PI);      //inc
      pm[i++] = 40.0;      //opn
      pm[i++] = 0.5;       //k
      pm[i++] = 2.0;       //gam
      pm[i++] = 0.5;      //xi
      pm[i++] = log(2.0); //mbh
      pm[i++] = 0.5;      //fellip
      pm[i++] = 0.5;      //fflow
      pm[i++] = log(0.01);
      pm[i++] = log(0.1);
      pm[i++] = log(0.01);
      pm[i++] = log(0.1);
      pm[i++] = 0.0;      //theta_rot
      pm[i++] = log(0.001); //sig_turb
      break;

    case 6:
      i=0;
      pm[i++] = log(4.0);   // mu
      pm[i++] = 1.0;        // beta
      pm[i++] = 0.25;        // F
      pm[i++] = cos(20.0/180.0*PI); // inc
      pm[i++] = 40.0;       // opn
      pm[i++] = 0.0;        // kappa
      pm[i++] = 1.0;        // gamma
      pm[i++] = 1.0;        // obscuration
      pm[i++] = log(2.0);  //mbh
      pm[i++] = 0.5;       //fellip
      pm[i++] = 0.4;       //fflow
      pm[i++] = log(0.01); //
      pm[i++] = log(0.1);  //
      pm[i++] = log(0.01); //
      pm[i++] = log(0.1);  //
      pm[i++] = 0.0;       // theta_rot
      pm[i++] = log(0.001);  // sig_turb
      break;

    case 7:
      i=0;
      pm[i++] = log(4.0);  //mu
      pm[i++] = 0.8;       //beta
      pm[i++] = 0.1;       //F
      pm[i++] = cos(20.0/180.0*PI);//inc
      pm[i++] = 40.0;      //opn
      pm[i++] = 0.0;       //kappa
      pm[i++] = 1.0;       //gamma
      pm[i++] = 0.0;       //xi

      pm[i++] = 0.5;      //fsh
      pm[i++] = log(8.0); //mu_un
      pm[i++] = 1.2;      //beta_un
      pm[i++] = 0.1;      //F_un
      pm[i++] = 20.0;     //opn_un

      pm[i++] = log(2.0); //mbh
      pm[i++] = 0.5;      //fellip
      pm[i++] = 0.4;      //fflow
      pm[i++] = log(0.01); 
      pm[i++] = log(0.1);
      pm[i++] = log(0.01);
      pm[i++] = log(0.1);
      pm[i++] = 0.0;      //theta_rot
      pm[i++] = 0.5;      //fellip_un
      pm[i++] = 0.4;      //fflow_un
      pm[i++] = log(0.001); // sig_turb
      break;

    case 8:
      i=0;
      pm[i++] = 50.0;   //theta_min
      pm[i++] = 20.0;   //dtheta_max
      pm[i++] = log(1.0); // r_min
      pm[i++] = log(4.0);      // fr_max
      pm[i++] = 1.0;      // gamma
      pm[i++] = 1.0;      // alpha
      pm[i++] = -2.0;      // lambda
      pm[i++] = -0.5;      // k
      pm[i++] = 0.0;      // xi
      pm[i++] = log(30.0);     // Rv
      pm[i++] = log(20.0);     // Rblr, should larger than r_max=r_max * fr_max
      pm[i++] = cos(30.0/180.0*PI);     // inc
      pm[i++] = log(4.0); // mbh
      break;
    
    case 9:
      i=0;
      pm[i++] = log(4.0);   //mu
      pm[i++] = 1.0;        //beta
      pm[i++] = 0.2;        //F
      pm[i++] = cos(30.0/180.0*PI);     // inc
      pm[i++] = 30.0;       //opn
      pm[i++] = log(4.0);   // mbh
      break;
  }
  return;
}

/*
 * print out the parameter values
 */
void print_par_value_sim(double *pm, int np)
{
  if(BLRmodel_name == NULL)
    return;
  
  int i;
  for(i=0; i<np; i++)
  {
    printf("%02d %15s = %f\n", i, BLRmodel_name[i], pm[i]);
  }
  return;
}

/*!
 * read in parameter values from a file fname
 */
void read_param_value(double *pm, char *fname)
{
  FILE *fp;
  int i;

  fp = fopen(fname, "r");
  if(fp == NULL)
  {
    printf("Error: cannot open file %s!\n", fname);
    exit(-1);
  }
  
  printf("%d\n", num_params_blr_model);

  for(i=0; i<num_params_blr_model; i++)
  {
    if(fscanf(fp, "%lf\n", &pm[i])<1)
    {
      printf("Error: unable to read in %d-th value in %s!\n", i, fname);
      exit(-1);
    }
    printf("%f\n", pm[i]);
  }

  fclose(fp);
  return;
}