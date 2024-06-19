/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifdef SpecAstro

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <fftw3.h>

#include "brains.h"

/*
 * calculate SA phase and 2D TF
 */
void calculate_sa_transfun2d_from_blrmodel(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i;
  double *sa_model, *rm_model;

  /* share all parameters */
  if(parset.flag_sa_par_mutual == 0)
  {
    gen_cloud_sample(pm, 5, 0);
    transfun_2d_cal_with_sample(transv, trans2d, n_vel);
    calculate_sa_with_sample(pm);
  }
  /* share mass and inclination */
  else if(parset.flag_sa_par_mutual == 1)
  {
    rm_model = (double *)pm;
    sa_model = rm_model + num_params_blr;
    sa_model[idx_sa_par_mutual[0]] = rm_model[idx_rm_par_mutual[0]]; //mbh
    sa_model[idx_sa_par_mutual[1]] = rm_model[idx_rm_par_mutual[1]]; //inc

    gen_cloud_sample(pm, 2, flag_save);
    transfun_2d_cal_with_sample(transv, trans2d, n_vel);

    gen_sa_cloud_sample((void *)sa_model, 3, flag_save);
    calculate_sa_with_sample(pm);
  }
  /* share inclination, opening angle and all dynamical parameters */
  else 
  {
    rm_model = (double *)pm;
    sa_model = rm_model + num_params_blr;

    sa_model[idx_sa_par_mutual[1]] = rm_model[idx_rm_par_mutual[1]]; //inc
    sa_model[idx_sa_par_mutual[2]] = rm_model[idx_rm_par_mutual[2]]; //opening angle
    /* dynamical parameters from mass on */
    for(i=idx_sa_par_mutual[0]; i<num_params_sa_blr_model; i++)
    {
      sa_model[i] = rm_model[i]; /* note the BLR models are identical */
    }

    gen_cloud_sample(pm, 2, flag_save);
    transfun_2d_cal_with_sample(transv, trans2d, n_vel);

    gen_sa_cloud_sample((void *)sa_model, 3, flag_save);
    calculate_sa_with_sample(pm);
  }
}

/*
 * calculate SA phase and 1D TF
 */
void calculate_sa_transfun_from_blrmodel(const void *pm, int flag_save)
{
  double *sa_model, *rm_model;

  if(parset.flag_sa_par_mutual == 0)
  {
    gen_cloud_sample(pm, 4, flag_save);
    transfun_1d_cal_with_sample();
    calculate_sa_with_sample(pm);
  }
  else
  {
    rm_model = (double *)pm;
    sa_model = rm_model + num_params_blr;
    sa_model[idx_sa_par_mutual[1]] = rm_model[idx_rm_par_mutual[1]]; // inc
    sa_model[idx_sa_par_mutual[2]] = rm_model[idx_rm_par_mutual[2]]; // opening angle

    gen_cloud_sample(pm, 1, flag_save);
    transfun_1d_cal_with_sample();

    gen_sa_cloud_sample((void *)sa_model, 3, flag_save);
    calculate_sa_with_sample(pm);
  }
}

/* 
 * calculate SA phase and line profile from cloud sample
 */
void calculate_sa_with_sample(const void *pm)
{
  int i, j, k, idV;
  double V, V_offset, dV, y, z, alpha, beta, flux_norm, phase, *phase_norm, *alpha_cent, *beta_cent;
  double DA, PA, FA, CO, cos_PA, sin_PA, alphac, betac;
  double *pmodel = (double *)pm;

  phase_norm = workspace_phase;
  alpha_cent = phase_norm + n_vel_sa_data;
  beta_cent = alpha_cent + n_vel_sa_data;

  /* angular size distance */
  DA = exp(pmodel[num_params_blr + num_params_sa_blr_model]);   
  /* position angle */         
  PA = pmodel[num_params_blr + num_params_sa_blr_model + 1]/180.0 * PI;  
  /* line flux scaling factor */
  FA = exp(pmodel[num_params_blr + num_params_sa_blr_model+2]);   
  /* 
   * line center offset: V+dV = (w-(w0+dw0))/(w0+dw0) ==> dV = - dw0/w0 
   * equivalent to redshift offset:  dz = -(1+z) dw0/w0
   */       
  CO = -pmodel[num_params_blr + num_params_sa_blr_model+3]/parset.sa_linecenter * C_Unit;  
  /* continuum offset */
  alphac = pmodel[num_params_blr + num_params_sa_blr_model + 4];
  betac  = pmodel[num_params_blr + num_params_sa_blr_model + 5];

  /* North through East (E of N), 180 deg ambiguity */
  cos_PA = cos(PA);
  sin_PA = sin(PA);

  dV = vel_sa_data[1] - vel_sa_data[0];
  
  for(i=0; i<n_vel_sa_data; i++)
  {
    for(k=0; k<n_base_sa_data; k++)
    {
      phase_sa[k*n_vel_sa_data + i] = 0.0;
    }
    phase_norm[i] = 0.0;
    alpha_cent[i] = 0.0;
    beta_cent[i] = 0.0;
  }

  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    y = clouds_alpha[i];
    z = clouds_beta[i];
    
    //alpha = y * cos_PA + z * sin_PA;
    //beta = -y * sin_PA + z * cos_PA;
    alpha = y;
    beta = z;
    
    for(j=0; j< parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i*parset.n_vel_per_cloud + j] + CO;
      V_offset = V + bin_offset * dV; /* bin type: center or left edge */
      if(V_offset<vel_sa_data[0] || V_offset >= vel_sa_data[n_vel_sa_data-1]+dV)
        continue;
      idV = (V_offset - vel_sa_data[0])/dV;
      
      phase_norm[idV] += clouds_weight[i];
      alpha_cent[idV] += alpha * clouds_weight[i];
      beta_cent[idV] += beta * clouds_weight[i];
    }
  }
  
  /* normalize line spectrum */
  flux_norm = 0.0;
  for(j=0; j<n_vel_sa_data; j++)
  {
    flux_norm += phase_norm[j];
  }
  flux_norm /= (sa_flux_norm * n_vel_sa_data);
  for(j=0; j<n_vel_sa_data; j++)
  {
    Fline_sa[j] = FA * phase_norm[j]/(flux_norm + EPS);
  }

  for(j=0; j<n_vel_sa_data; j++)
  {
    y = alpha_cent[j];
    z = beta_cent[j];

    alpha_cent[j] = y * cos_PA + z * sin_PA;
    beta_cent[j] = -y * sin_PA + z * cos_PA;

    alpha_cent[j] = (alpha_cent[j]/(phase_norm[j]+EPS) - alphac) / DA;
    beta_cent[j]  =  (beta_cent[j]/(phase_norm[j]+EPS) - betac) / DA;
  }
  
  /* phi = -2*pi * f_line * B/lambda * X/DA for phase
   *     =         f_line * j * X/DA for photocenter
   *
   * note that phase and photocenter are scaled 
   */
  for(k=0; k<n_base_sa_data; k++)
  {
    for(j=0; j<n_vel_sa_data; j++)
    { 
      phase = base_sa_data[k*2] * alpha_cent[j] + base_sa_data[k*2 + 1] * beta_cent[j];
      phase_sa[k*n_vel_sa_data + j] = sign * Fline_sa[j]/(1.0+Fline_sa[j]) * phase;
    }
  }

  /* smooth */
  sa_smooth_run(vel_sa_data, Fline_sa, n_vel_sa_data, phase_sa, n_base_sa_data);

  return;
}

/* 
 * calculate SA phase and line profile from cloud sample at simulated mesh.
 */
void calculate_sa_sim_with_sample(const void *pm, double *vel_sa, int n_vel_sa, double *base_sa, int n_base_sa,
                                  double *p_sa, double *F_sa)
{
  int i, j, k, idV;
  double V, V_offset, dV, y, z, alpha, beta, flux_norm, phase, *phase_norm, *alpha_cent, *beta_cent;
  double DA, PA, FA, CO, cos_PA, sin_PA, alphac, betac;
  double *pmodel = (double *)pm;

  phase_norm = workspace_phase;
  alpha_cent = phase_norm + n_vel_sa;
  beta_cent = alpha_cent + n_vel_sa;

  /* angular size distance */
  DA = exp(pmodel[num_params_blr + num_params_sa_blr_model]);   
  /* position angle */         
  PA = pmodel[num_params_blr + num_params_sa_blr_model + 1]/180.0 * PI;  
  /* line flux scaling factor */
  FA = exp(pmodel[num_params_blr + num_params_sa_blr_model+2]);   
  /* 
   * line center offset: V+dV = (w-(w0+dw0))/(w0+dw0) ==> dV = - dw0/w0 
   * equivalent to redshift offset:  dz = -(1+z) dw0/w0
   */       
  CO = -pmodel[num_params_blr + num_params_sa_blr_model+3]/parset.sa_linecenter * C_Unit; 
  /* continuum offset */
  alphac = pmodel[num_params_blr + num_params_sa_blr_model + 4];
  betac  = pmodel[num_params_blr + num_params_sa_blr_model + 5];

  /* North through East (E of N), 180 deg ambiguity */
  cos_PA = cos(PA);
  sin_PA = sin(PA);

  dV = vel_sa[1] - vel_sa[0];
  
  for(i=0; i<n_vel_sa; i++)
  {
    for(k=0; k<n_base_sa; k++)
    {
      p_sa[k*n_vel_sa + i] = 0.0;
    }
    phase_norm[i] = 0.0;
    alpha_cent[i] = 0.0;
    beta_cent[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    y = clouds_alpha[i];
    z = clouds_beta[i];
    
    //alpha = y * cos_PA + z * sin_PA;
    //beta = -y * sin_PA + z * cos_PA;
    alpha = y;
    beta = z;
    
    for(j=0; j< parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i*parset.n_vel_per_cloud + j] + CO;
      V_offset = V + bin_offset * dV; /* bin type: center or left edge */
      if(V_offset<vel_sa[0] || V_offset >= vel_sa[n_vel_sa-1]+dV)
        continue;
      idV = (V_offset - vel_sa[0])/dV;
      
      phase_norm[idV] += clouds_weight[i];
      alpha_cent[idV] += alpha * clouds_weight[i];
      beta_cent[idV] += beta * clouds_weight[i];
    }
  }

  /* normalize line spectrum */
  flux_norm = 0.0;
  for(j=0; j<n_vel_sa; j++)
  {
    flux_norm += phase_norm[j];
  }
  flux_norm /= (sa_flux_norm * n_vel_sa);
  for(j=0; j<n_vel_sa; j++)
  {
    F_sa[j] = FA * phase_norm[j]/(flux_norm + EPS);
  }

  for(j=0; j<n_vel_sa; j++)
  {
    y = alpha_cent[j];
    z = beta_cent[j];

    alpha_cent[j] = y * cos_PA + z * sin_PA;
    beta_cent[j] = -y * sin_PA + z * cos_PA;

    alpha_cent[j] = (alpha_cent[j]/(phase_norm[j]+EPS) - alphac) / DA;
    beta_cent[j]  = (beta_cent[j]/(phase_norm[j]+EPS)  - betac) / DA;
  }
  
  /* phi = -2*pi * f_line * B/lambda * X/DA for phase
   *     =         f_line * j * X/DA for photocenter
   *
   * note that phase and photocenter are scaled 
   */
  for(k=0; k<n_base_sa; k++)
  {
    for(j=0; j<n_vel_sa; j++)
    { 
      phase = base_sa[k*2] * alpha_cent[j] + base_sa[k*2 + 1] * beta_cent[j];
      p_sa[k*n_vel_sa + j] = sign * F_sa[j]/(1.0+F_sa[j]) * phase;
    }
  }

  /* smooth */
  sa_smooth_run(vel_sa, F_sa, n_vel_sa, p_sa, n_base_sa);
  
  return;
}

/* 
 * calculate SA phase and line profile from a BLR model.
 */
void calculate_sa_from_blrmodel(const void *pm, int flag_save)
{
  gen_sa_cloud_sample(pm, 3, flag_save);
  calculate_sa_with_sample(pm);
  return;
}

void set_sa_blr_model()
{
  switch(parset.flag_sa_blrmodel)
  {
    case 0:
      num_params_sa_blr_model = num_params_MyBLRmodel2d;
      gen_sa_cloud_sample = gen_cloud_sample_mymodel;
      BLRmodel_sa_name = MyBLRmodel_name;
      break;

    case 1:
      num_params_sa_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model1;
      BLRmodel_sa_name = BLRmodel1_name;
      break;

    case 2:
      num_params_sa_blr_model = sizeof(BLRmodel2)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model2;
      BLRmodel_sa_name = BLRmodel2_name;
      break;

    case 3:
      num_params_sa_blr_model = sizeof(BLRmodel3)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model3;
      BLRmodel_sa_name = BLRmodel3_name;
      break;

    case 4:
      num_params_sa_blr_model = sizeof(BLRmodel4)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model4;
      BLRmodel_sa_name = BLRmodel4_name;
      break;

    case 5:
      num_params_sa_blr_model = sizeof(BLRmodel5)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model5;
      BLRmodel_sa_name = BLRmodel5_name;
      break;
    
    case 6:
      num_params_sa_blr_model = sizeof(BLRmodel6)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model6;
      BLRmodel_sa_name = BLRmodel6_name;
      break;
    
    case 7:
      num_params_sa_blr_model = sizeof(BLRmodel7)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model7;
      BLRmodel_sa_name = BLRmodel7_name;
      break;

    case 8:
      num_params_sa_blr_model = sizeof(BLRmodel8)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model8;
      BLRmodel_sa_name = BLRmodel8_name;
      break;
    
    case 9:
      num_params_sa_blr_model = sizeof(BLRmodel9)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model9;
      BLRmodel_sa_name = BLRmodel9_name;
      break;

    default:
      num_params_sa_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model1;
      BLRmodel_sa_name = BLRmodel1_name;
      break;
  }

  return;
}

/*
 * set the index of mutual parameters (BH mass and inclination) in RM and SA BLR models.
 */
void set_idx_par_mutual()
{
  switch(parset.flag_blrmodel)
  {
    case 0:
      idx_rm_par_mutual[0] = offsetof(MyBLRmodel, mbh);
      idx_rm_par_mutual[1] = offsetof(MyBLRmodel, inc);
      idx_rm_par_mutual[2] = offsetof(MyBLRmodel, opn);
      break;

    case 1:
      idx_rm_par_mutual[0] = offsetof(BLRmodel1, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel1, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel1, opn);
      break;

    case 2:
      idx_rm_par_mutual[0] = offsetof(BLRmodel2, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel2, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel2, opn);
      break;
    
    case 3:
      idx_rm_par_mutual[0] = offsetof(BLRmodel3, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel3, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel3, opn);
      break;

    case 4:
      idx_rm_par_mutual[0] = offsetof(BLRmodel4, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel4, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel4, opn);
      break;
    
    case 5:
      idx_rm_par_mutual[0] = offsetof(BLRmodel5, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel5, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel5, opn);
      break;

    case 6:
      idx_rm_par_mutual[0] = offsetof(BLRmodel6, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel6, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel6, opn);
      break;
    
    case 7:
      idx_rm_par_mutual[0] = offsetof(BLRmodel7, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel7, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel7, opn);
      break;
    
    case 8:
      idx_rm_par_mutual[0] = offsetof(BLRmodel8, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel8, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel8, theta_min);
      break;
    
    case 9:
      idx_rm_par_mutual[0] = offsetof(BLRmodel9, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel9, inc);
      idx_rm_par_mutual[2] = offsetof(BLRmodel9, opn);
      break;
  }

  switch(parset.flag_sa_blrmodel)
  {
    case 0:
      idx_sa_par_mutual[0] = offsetof(MyBLRmodel, mbh);
      idx_sa_par_mutual[1] = offsetof(MyBLRmodel, inc);
      idx_sa_par_mutual[2] = offsetof(MyBLRmodel, opn);
      break;

    case 1:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel1, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel1, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel1, opn);
      break;

    case 2:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel2, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel2, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel2, opn);
      break;
    
    case 3:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel3, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel3, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel3, opn);
      break;

    case 4:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel4, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel4, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel4, opn);
      break;
    
    case 5:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel5, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel5, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel5, opn);
      break;

    case 6:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel6, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel6, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel6, opn);
      break;
    
    case 7:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel7, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel7, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel7, opn);
      break;
    
    case 8:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel8, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel8, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel8, theta_min);
      break;
    
    case 9:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel9, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel9, inc);
      idx_sa_par_mutual[2] = offsetof(SABLRmodel9, opn);
      break;
  }

  idx_sa_par_mutual[0] /= sizeof(double);
  idx_sa_par_mutual[1] /= sizeof(double);
  idx_sa_par_mutual[2] /= sizeof(double);
  idx_rm_par_mutual[0] /= sizeof(double);
  idx_rm_par_mutual[1] /= sizeof(double);
  idx_rm_par_mutual[2] /= sizeof(double);

  return;
}

/*
 * setup  SA BLR model parameter range. 
 */
// model 1
void set_sa_blr_range_model1()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //lambda
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.5;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //sa_blr_range_model[2][1] = fmin(sa_blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 2
void set_sa_blr_range_model2()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //sigr
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigtheta
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  return;
}

// model 3
void set_sa_blr_range_model3()
{
  int i;
  
  i = 0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn 
  sa_blr_range_model[i][0] = 0.0;   // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[3][1] + sa_blr_range_model[4][1]));
  return;
}

// model 4
void set_sa_blr_range_model4()
{
  int i;
  
  i = 0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = sqrt(2.0)/2.0;
  //q
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[3][1] + sa_blr_range_model[4][1]));
  return;
}

// model 5
void set_sa_blr_range_model5()
{
  int i;
  
  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //Fin
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //Fout
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(10.0);
  //alpha
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 3.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //beta
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;
  //sig_turb
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);

  //rcloud_max_set = fmax(rcloud_max_set, exp(sa_blr_range_model[2][1] + sa_blr_range_model[4][1]));

  return;
}


// model 6
void set_sa_blr_range_model6()
{
  int i;
  
  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //inc
  sa_blr_range_model[i][0] = 0.0; //in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //gamma
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;
  //sig_turb
  sa_blr_range_model[i][0] = log(0.0001);
  sa_blr_range_model[i++][1] = log(0.1);

  //rcloud_max_set = 44.85;
  //sa_blr_range_model[2][1] = fmin(sa_blr_range_model[2][1], log(rcloud_max_set));
  return;
}

// model 7
void set_sa_blr_range_model7()
{
  int i;
  

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in rad
  sa_blr_range_model[i++][1] = 90.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //gamma
  sa_blr_range_model[i][0] = 1.0;
  sa_blr_range_model[i++][1] = 5.0;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //fsh
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //mu_un
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta_un
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //opn_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;

  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);
  //fellip
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //sigr_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_circ
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //sigr_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);
  //sigthe_rad
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(1.0);
  //theta_rot
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 90.0;

  //fellip_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //fflow_un
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;

  //sig_turb
  sa_blr_range_model[i][0] = log(0.001);
  sa_blr_range_model[i++][1] = log(0.1);

  return;
}

// model 8
void set_sa_blr_range_model8()
{
  int i;
  
  i = 0;
  //theta_min
  sa_blr_range_model[i][0] = 20.0;     // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //dtheta_max
  sa_blr_range_model[i][0] = 0.0;      // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //r_min
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.25);
  //fr_max
  sa_blr_range_model[i][0] = 1.0;  //
  sa_blr_range_model[i++][1] = 10.0;
  //gamma
  sa_blr_range_model[i][0] = 0.0;  //
  sa_blr_range_model[i++][1] = 3.0;
  //alpha
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 3.0;
  //lambda
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 0.0;
  //k
  sa_blr_range_model[i][0] = -0.5;
  sa_blr_range_model[i++][1] = 0.5;
  //xi
  sa_blr_range_model[i][0] = 0.0;
  sa_blr_range_model[i++][1] = 1.0;
  //Rv
  sa_blr_range_model[i][0] = log(10.0);
  sa_blr_range_model[i++][1] = log(50.0);
  //Rblr
  sa_blr_range_model[i][0] = log(10.0);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //inc
  sa_blr_range_model[i][0] = 0.0;   // cos(inc)
  sa_blr_range_model[i++][1] = 1.0;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);

  return;
}

// model 9
void set_sa_blr_range_model9()
{
  int i;

  i = 0;
  //mu
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //beta
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 2.0;
  //F
  sa_blr_range_model[i][0] = 0.001;
  sa_blr_range_model[i++][1] = 0.999;
  //inc
  sa_blr_range_model[i][0] = 0.0;  // in cosine
  sa_blr_range_model[i++][1] = 1.0;
  //opn
  sa_blr_range_model[i][0] = 0.0;  // in degree
  sa_blr_range_model[i++][1] = 90.0;
  //mbh
  sa_blr_range_model[i][0] = log(mass_range[0]);
  sa_blr_range_model[i++][1] = log(mass_range[1]);

  return;
}

/*!
 * This function copes with parameter fixing.\n
 * Only fix BLR model parameters. 
 */
void set_par_fix_sa_blrmodel()
{
  int i;
  char *pstr;
  
  npar_fix = 0;

  if(thistask == roottask)
  {
    pstr = parset.sa_str_par_fix_val;
    // set the default value if not provided.
    for(i=strlen(parset.sa_str_par_fix); i<num_params_sa_blr_model; i++)
      parset.sa_str_par_fix[i] = '0';

    for(i=0; i<num_params_sa_blr_model; i++)
    {
      if(parset.sa_str_par_fix[i] == '0')
      {
        par_fix[num_params_blr + i] = 0;
        par_fix_val[num_params_blr + i] = -DBL_MAX;  /* set to be the smallest value */
      }
      else if(parset.sa_str_par_fix[i] == '1')
      {
        if(pstr == NULL)
        {
          printf("# %d-th SA BLR parameter value is not provided (counting from 0).\n", i);
          exit(0);
        }
        par_fix[num_params_blr + i] = 1;
        sscanf(pstr, "%lf", &par_fix_val[num_params_blr + i]);
        npar_fix++;
        printf("# %d-th SA BLR parameter fixed, value= %f.\n", i, par_fix_val[num_params_blr + i]);
        pstr = strchr(pstr, ':'); /* values are separated by ":" */
        if(pstr!=NULL)
        {
          pstr++;
        }
      }
      else   // default value
      {
        par_fix[num_params_blr + i] = 0;
        par_fix_val[num_params_blr + i] = -DBL_MAX;
      }
    }
  }

  MPI_Bcast(par_fix + num_params_blr, num_params_sa_blr_model, MPI_INT, roottask, MPI_COMM_WORLD);
  MPI_Bcast(par_fix_val + num_params_blr, num_params_sa_blr_model, MPI_DOUBLE, roottask, MPI_COMM_WORLD);

  return;
}


int nd_fft_sa, nd_fft_sa_cal, npad_sa;

fftw_complex *sa_data_fft, *sa_resp_fft, *sa_conv_fft;
fftw_plan sa_pdata, sa_presp, sa_pback;
double *sa_real_data, *sa_real_resp, *sa_real_conv;

/*!
 * This function initiates workspace for FFT.
 */
void sa_smooth_init(int n_v_sa, const double *v_sa, double sigV)
{
  int i;
  double dV;

  npad_sa = fmin(n_v_sa * 0.1, 20);
  nd_fft_sa = n_v_sa+npad_sa*2;
  nd_fft_sa_cal = nd_fft_sa/2+1;

  sa_data_fft = (fftw_complex *) fftw_malloc((nd_fft_sa_cal) * sizeof(fftw_complex));
  sa_resp_fft = (fftw_complex *) fftw_malloc((nd_fft_sa_cal) * sizeof(fftw_complex));
  sa_conv_fft = (fftw_complex *) fftw_malloc((nd_fft_sa_cal) * sizeof(fftw_complex));

  sa_real_data = (double *)fftw_malloc(nd_fft_sa * sizeof(double));
  //real_resp = (double *)fftw_malloc(nd_fft_sa * sizeof(double));
  sa_real_conv = (double *)fftw_malloc(nd_fft_sa * sizeof(double));

  sa_pdata = fftw_plan_dft_r2c_1d(nd_fft_sa, sa_real_data, sa_data_fft, FFTW_PATIENT);
  //presp = fftw_plan_dft_r2c_1d(nd_fft_sa, sa_real_resp, sa_resp_fft, FFTW_PATIENT);
  sa_pback = fftw_plan_dft_c2r_1d(nd_fft_sa, sa_conv_fft, sa_real_conv, FFTW_PATIENT);

  dV = v_sa[1] - v_sa[0];
  for(i=0; i<nd_fft_sa_cal; i++)
  {
    sa_resp_fft[i][0] = exp(-2.0 * PI*PI * sigV/dV*sigV/dV * i*i*1.0/nd_fft_sa/nd_fft_sa)/nd_fft_sa;
    sa_resp_fft[i][1] = 0.0;
  }
  return;
}

/*!
 * This function finalizes FFT.
 */
void sa_smooth_end()
{
  fftw_destroy_plan(sa_pdata);
  fftw_destroy_plan(sa_pback);

  fftw_free(sa_data_fft);
  fftw_free(sa_resp_fft);
  fftw_free(sa_conv_fft);
  
  fftw_free(sa_real_data);
  //fftw_free(sa_real_resp);
  fftw_free(sa_real_conv);
  return;
}

void sa_smooth_run(double *v_sa, double *F_sa, int n_v_sa, double *p_sa, int n_base_sa)
{
  int i, j;
  /* first profile */
  memcpy(sa_real_data+npad_sa, F_sa, n_v_sa*sizeof(double));
  for(i=0; i<npad_sa; i++)
    sa_real_data[i] = sa_real_data[nd_fft_sa-1-i] = 0.0;
  
  /* FFT of line */
  fftw_execute(sa_pdata);
  /* complex multiply and inverse FFT 
   * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
   * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
   */
  for(i=0; i<nd_fft_sa_cal; i++)
  {
    sa_conv_fft[i][0] = sa_data_fft[i][0]*sa_resp_fft[i][0] - sa_data_fft[i][1]*sa_resp_fft[i][1];
    sa_conv_fft[i][1] = sa_data_fft[i][0]*sa_resp_fft[i][1] + sa_data_fft[i][1]*sa_resp_fft[i][0];
  }
  fftw_execute(sa_pback);
  memcpy(F_sa, sa_real_conv+npad_sa, n_v_sa*sizeof(double));

  /* then phase */
  for(j=0; j<n_base_sa; j++)
  {
    memcpy(sa_real_data+npad_sa, &p_sa[j*n_v_sa], n_v_sa*sizeof(double));
    for(i=0; i<npad_sa; i++)
      sa_real_data[i] = sa_real_data[nd_fft_sa-1-i] = 0.0;
    
    /* FFT of line */
    fftw_execute(sa_pdata);
    /* complex multiply and inverse FFT 
     * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
     * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
     */
    for(i=0; i<nd_fft_sa_cal; i++)
    {
      sa_conv_fft[i][0] = sa_data_fft[i][0]*sa_resp_fft[i][0] - sa_data_fft[i][1]*sa_resp_fft[i][1];
      sa_conv_fft[i][1] = sa_data_fft[i][0]*sa_resp_fft[i][1] + sa_data_fft[i][1]*sa_resp_fft[i][0];
    }
    fftw_execute(sa_pback);
    memcpy(&p_sa[j*n_v_sa], sa_real_conv+npad_sa, n_v_sa*sizeof(double));
  }
  return;
}
#endif