/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifdef SA

#include <math.h>
#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "brains.h"

void calculate_sa_transfun2d_from_blrmodel(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  double *sa_model, *rm_model;

  if(parset.flag_sa_par_mutual == 0)
  {
    gen_cloud_sample(pm, 5, 0);
    transfun_2d_cal_with_sample(transv, trans2d, n_vel);
    calculate_sa_with_sample(pm);
  }
  else
  {
    rm_model = (double *)pm;
    sa_model = rm_model + num_params_blr;
    sa_model[idx_sa_par_mutual[0]] = rm_model[idx_rm_par_mutual[0]];
    sa_model[idx_sa_par_mutual[1]] = rm_model[idx_rm_par_mutual[1]];

    gen_cloud_sample(pm, 2, flag_save);
    transfun_2d_cal_with_sample(transv, trans2d, n_vel);

    gen_sa_cloud_sample((void *)sa_model, 3, flag_save);
    calculate_sa_with_sample(pm);
  }
}

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
    sa_model[idx_sa_par_mutual[0]] = rm_model[idx_rm_par_mutual[0]];
    sa_model[idx_sa_par_mutual[1]] = rm_model[idx_rm_par_mutual[1]];

    gen_cloud_sample(pm, 1, flag_save);
    transfun_1d_cal_with_sample();

    gen_sa_cloud_sample((void *)sa_model, 3, flag_save);
    calculate_sa_with_sample(pm);
  }
}

/* 
 * calculate SA phase and line profile.
 */
void calculate_sa_with_sample(const void *pm)
{
  int i, j, k, idV;
  double V, dV, y, z, alpha, beta, flux_norm, phase, *phase_norm, *alpha_cent, *beta_cent;
  double DA, PA, FA, CO, cos_PA, sin_PA, Vres;
  double *pmodel = (double *)pm;

  phase_norm = workspace_phase;
  alpha_cent = phase_norm + n_vel_sa_data;
  beta_cent = alpha_cent + n_vel_sa_data;

  /* angular size distance */
  DA = exp(pmodel[num_params_blr + num_params_sa_blr_model]);   
  /* position angle */         
  PA = pmodel[num_params_blr + num_params_sa_blr_model + 1]/180.0 * PI;  
  /* line flux scaling center */
  FA = exp(pmodel[num_params_blr + num_params_sa_blr_model+2]);   
  /* 
   * line center offset: V+dV = (w-(w0+dw0))/(w0+dw0) ==> dV = - dw0/w0 
   * equivalent to redshift offset:  dz = -(1+z) dw0/w0
   */       
  CO = -pmodel[num_params_blr + num_params_sa_blr_model+3]/parset.sa_linecenter * C_Unit;  

  /* instrument broadening */
  Vres = parset.sa_InstRes; 

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
    
    alpha = y * cos_PA + z * sin_PA;
    beta = -y * sin_PA + z * cos_PA;

    for(j=0; j< parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i*parset.n_vel_per_cloud + j] + CO + Vres * gsl_ran_ugaussian(gsl_r);
      if(V<vel_sa_data[0] || V >= vel_sa_data[n_vel_sa_data-1]+dV)
        continue;
      idV = (V - vel_sa_data[0])/dV;
      
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
    alpha_cent[j] = (alpha_cent[j]/(phase_norm[j]+EPS)) / DA;
    beta_cent[j] =  (beta_cent[j]/(phase_norm[j]+EPS)) / DA;
  }
  
  /* phi = -2*pi * f_line * B/lambda * X/DA */
  for(k=0; k<n_base_sa_data; k++)
  {
    for(j=0; j<n_vel_sa_data; j++)
    { 
      phase = base_sa_data[k*2] * alpha_cent[j] + base_sa_data[k*2 + 1] * beta_cent[j];
      phase_sa[k*n_vel_sa_data + j] = -Fline_sa[j]/(1.0+Fline_sa[j]) * phase;
    }
  }

  return;
}

/* 
 * calculate SA phase and line profile.
 */
void calculate_sa_from_blrmodel(const void *pm, int flag_save)
{
  int i, j, k, idV;
  double V, dV, y, z, alpha, beta, flux_norm, phase, *phase_norm, *alpha_cent, *beta_cent;
  double DA, PA, FA, CO, cos_PA, sin_PA, Vres;
  double *pmodel = (double *)pm;

  phase_norm = workspace_phase;
  alpha_cent = phase_norm + n_vel_sa_data;
  beta_cent = alpha_cent + n_vel_sa_data;

  /* angular size distance */
  DA = exp(pmodel[num_params_blr + num_params_sa_blr_model]);   
  /* position angle */         
  PA = pmodel[num_params_blr + num_params_sa_blr_model + 1]/180.0 * PI;  
  /* line flux scaling center */
  FA = exp(pmodel[num_params_blr + num_params_sa_blr_model+2]);   
  /* 
   * line center offset: V+dV = (w-(w0+dw0))/(w0+dw0) ==> dV = - dw0/w0 
   * equivalent to redshift offset:  dz = -(1+z) dw0/w0
   */       
  CO = -pmodel[num_params_blr + num_params_sa_blr_model+3]/parset.sa_linecenter * C_Unit;  

  /* instrument broadening */
  Vres = parset.sa_InstRes; 

  cos_PA = cos(PA);
  sin_PA = sin(PA);

  gen_sa_cloud_sample(pm, 3, 0);

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
    
    alpha = y * cos_PA + z * sin_PA;
    beta = -y * sin_PA + z * cos_PA;

    for(j=0; j< parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i*parset.n_vel_per_cloud + j] + CO + Vres * gsl_ran_ugaussian(gsl_r);
      if(V<vel_sa_data[0] || V >= vel_sa_data[n_vel_sa_data-1]+dV)
        continue;
      idV = (V - vel_sa_data[0])/dV;
      
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
    alpha_cent[j] = (alpha_cent[j]/(phase_norm[j]+EPS)) / DA;
    beta_cent[j] =  (beta_cent[j]/(phase_norm[j]+EPS)) / DA;
  }
  
  /* phi = -2*pi * f_line * B/lambda * X/DA */
  for(k=0; k<n_base_sa_data; k++)
  {
    for(j=0; j<n_vel_sa_data; j++)
    { 
      phase = base_sa_data[k*2] * alpha_cent[j] + base_sa_data[k*2 + 1] * beta_cent[j];
      phase_sa[k*n_vel_sa_data + j] = -Fline_sa[j]/(1.0+Fline_sa[j]) * phase;
    }
  }

  return;
}

void set_sa_blr_model()
{
  switch(parset.flag_sa_blrmodel)
  {
    case 0:
      num_params_sa_blr_model = num_params_MyBLRmodel2d;
      gen_sa_cloud_sample = gen_cloud_sample_mymodel;
      break;

    case 1:
      num_params_sa_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model1;
      break;

    case 2:
      num_params_sa_blr_model = sizeof(BLRmodel2)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model2;
      break;

    case 3:
      num_params_sa_blr_model = sizeof(BLRmodel3)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model3;
      break;

    case 4:
      num_params_sa_blr_model = sizeof(BLRmodel4)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model4;
      break;

    case 5:
      num_params_sa_blr_model = sizeof(BLRmodel5)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model5;
      break;
    
    case 6:
      num_params_sa_blr_model = sizeof(BLRmodel6)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model6;
      break;
    
    case 7:
      num_params_sa_blr_model = sizeof(BLRmodel7)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model7;
      break;

    case 8:
      num_params_sa_blr_model = sizeof(BLRmodel8)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model8;
      break;
    
    case 9:
      num_params_sa_blr_model = sizeof(BLRmodel9)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model9;
      break;

    default:
      num_params_sa_blr_model = sizeof(BLRmodel1)/sizeof(double);
      gen_sa_cloud_sample = gen_cloud_sample_model1;
      break;
  }

  return;
}

void set_idx_par_mutual()
{
  switch(parset.flag_blrmodel)
  {
    case 0:
      idx_rm_par_mutual[0] = offsetof(MyBLRmodel, mbh);
      idx_rm_par_mutual[1] = offsetof(MyBLRmodel, inc);
      break;

    case 1:
      idx_rm_par_mutual[0] = offsetof(BLRmodel1, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel1, inc);
      break;

    case 2:
      idx_rm_par_mutual[0] = offsetof(BLRmodel2, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel2, inc);
      break;
    
    case 3:
      idx_rm_par_mutual[0] = offsetof(BLRmodel3, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel3, inc);
      break;

    case 4:
      idx_rm_par_mutual[0] = offsetof(BLRmodel4, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel4, inc);
      break;
    
    case 5:
      idx_rm_par_mutual[0] = offsetof(BLRmodel5, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel5, inc);
      break;

    case 6:
      idx_rm_par_mutual[0] = offsetof(BLRmodel6, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel6, inc);
      break;
    
    case 7:
      idx_rm_par_mutual[0] = offsetof(BLRmodel7, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel7, inc);
      break;
    
    case 8:
      idx_rm_par_mutual[0] = offsetof(BLRmodel8, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel8, inc);
      break;
    
    case 9:
      idx_rm_par_mutual[0] = offsetof(BLRmodel9, mbh);
      idx_rm_par_mutual[1] = offsetof(BLRmodel9, inc);
      break;
  }

  switch(parset.flag_sa_blrmodel)
  {
    case 0:
      idx_sa_par_mutual[0] = offsetof(MyBLRmodel, mbh);
      idx_sa_par_mutual[1] = offsetof(MyBLRmodel, inc);
      break;

    case 1:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel1, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel1, inc);
      break;

    case 2:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel2, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel2, inc);
      break;
    
    case 3:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel3, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel3, inc);
      break;

    case 4:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel4, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel4, inc);
      break;
    
    case 5:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel5, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel5, inc);
      break;

    case 6:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel6, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel6, inc);
      break;
    
    case 7:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel7, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel7, inc);
      break;
    
    case 8:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel8, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel8, inc);
      break;
    
    case 9:
      idx_sa_par_mutual[0] = offsetof(SABLRmodel9, mbh);
      idx_sa_par_mutual[1] = offsetof(SABLRmodel9, inc);
      break;
  }

  idx_sa_par_mutual[0] /= sizeof(double);
  idx_sa_par_mutual[1] /= sizeof(double);
  idx_rm_par_mutual[0] /= sizeof(double);
  idx_rm_par_mutual[1] /= sizeof(double);

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
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
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
  //alpha
  sa_blr_range_model[i][0] = -3.0;
  sa_blr_range_model[i++][1] = 3.0;
  //Rin
  sa_blr_range_model[i][0] = log(0.1);
  sa_blr_range_model[i++][1] = log(rcloud_max_set*0.5);
  //F
  sa_blr_range_model[i][0] = log(1.0);
  sa_blr_range_model[i++][1] = log(1.0e2);
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

#endif