/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*
 * this file is to calculate spectroastrometric reverberation mapping.
 */

#ifdef SpecAstro

#include <fftw3.h>

#include "brains.h"

/* 
 * calculate SA and spectral TF
 *
 */
void transfun_sarm_cal_cloud(const void *pm, double *transv, double *trans2d, double *trans_sarm_alpha, double *trans_sarm_beta, 
                             int n_vel, int flag_save)
{
  gen_cloud_sample(pm, 5, 0);
  transfun_sarm_cal_with_sample(transv, trans2d, trans_sarm_alpha, trans_sarm_beta, n_vel);
  return;
}

/*
 * calculate SA and spectral TF given the cloud sample
 *
 */
void transfun_sarm_cal_with_sample(double *transv, double *trans2d, double *trans_alpha, double *trans_beta, int n_vel)
{
  int i, j, idV, idt;
  double tau_min, tau_max;
  double V, dis, dV, dTransTau, V_offset, Anorm;

  tau_min = clouds_tau[0];
  tau_max = clouds_tau[0];
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    if(tau_min > clouds_tau[i])
      tau_min = clouds_tau[i];
    if(tau_max < clouds_tau[i])
      tau_max = clouds_tau[i];
  }

  dTransTau = (tau_max - tau_min)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = tau_min + dTransTau * i;
  }

  dV =(transv[1] - transv[0]); /* velocity grid width */
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0;j<n_vel;j++)
    {
      trans2d[i*n_vel+j]=0.0;   /* cleanup of transfer function */
      trans_alpha[i*n_vel+j] = 0.0;
      trans_beta[i*n_vel+j] = 0.0;
    }
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    dis = clouds_tau[i];
    idt = (dis - tau_min)/dTransTau;

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      V = clouds_vel[i*parset.n_vel_per_cloud + j];
      V_offset = V + bin_offset * dV; /* bin type: center or left edge */
      if(V_offset < transv[0] || V_offset >= transv[n_vel-1] + dV )
        continue;
      idV = (V_offset - transv[0])/dV; 
      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += clouds_weight[i];
      trans_alpha[idt*n_vel + idV] += clouds_weight[i] * clouds_alpha[i];
      trans_beta[idt*n_vel + idV] += clouds_weight[i] * clouds_beta[i];
    }
  }

  /* normalize transfer functions */
  Anorm = 0.0;
  for(i=0; i<parset.n_tau; i++)
    for(j=0; j<n_vel; j++)
    {
      Anorm += trans2d[i*n_vel+j];
    }
  Anorm *= (dV * dTransTau);
  
  Anorm += EPS;
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<n_vel; j++)
    {
      trans2d[i*n_vel+j] /= Anorm;
      trans_alpha[i*n_vel+j] /= Anorm;
      trans_beta[i*n_vel+j] /= Anorm;
    }
  }

  return;
}

/*
 * calculate SARM signal on simulated grids 
 * inputs:  pm, tline, vel_sa, trans2d, trans_alpha, trans_beta, n_vel, n_line,
 *          base, n_base
 * outputs: phase, fline, photocenter_alpha, photocenter_beta
 */
void calculate_sarm_sim_with_sample(const void *pm, double *tline, double *vel_sa, double *trans2d, 
                                    double *trans_alpha, double *trans_beta,
                                    int n_vel, int n_line, double *base, 
                                    int n_base, double *phase, double *fline, 
                                    double *momentum_alpha, double *momentum_beta,
                                    double *photocenter_alpha, double *photocenter_beta)
{
  int i, j, k, m;
  double tau, tl, tc, fcon_rm, dTransTau, ratio, flux_ratio, fcon;
  double DA, PA, cos_PA, sin_PA, y, z;
  double *pmodel = (double *)pm;

  dTransTau = TransTau[1] - TransTau[0];

  /* angular size distance */
  DA = exp(pmodel[num_params_blr + num_params_sa_blr_model]); 
  /* position angle */         
  PA = pmodel[num_params_blr + num_params_sa_blr_model + 1]/180.0 * PI;  

  /* North through East (E of N), 180 deg ambiguity */
  cos_PA = cos(PA);
  sin_PA = sin(PA);

  /* loop over time of line */
  for(j=0;j<n_line; j++)
  {
    for(i=0; i<n_vel; i++)
    {
      fline[j*n_vel + i] = 0.0;
      momentum_alpha[j*n_vel + i] = 0.0;
      momentum_beta[j*n_vel + i] = 0.0;
    }

    tl = tline[j];
    /* interpret to get the present continuum flux */
    fcon = Fcon_sarm[j];

    for(k=0; k<parset.n_tau; k++)
    {
      tau = TransTau[k];
      tc = tl - tau;
      fcon_rm = gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, tc, gsl_acc);

      for(i=0; i<n_vel; i++)
      {
        fline[j*n_vel + i] += trans2d[k*n_vel+i] * fcon_rm;
      
        momentum_alpha[j*n_vel + i] += trans_alpha[k*n_vel + i] * fcon_rm;
        momentum_beta[j*n_vel + i] += trans_beta[k*n_vel + i] * fcon_rm;
      }
    }

    for(i=0; i<n_vel; i++)
    {
      fline[j*n_vel + i] *= dTransTau;
      momentum_alpha[j*n_vel + i] *= dTransTau;
      momentum_beta[j*n_vel + i] *= dTransTau;

      y = momentum_alpha[j*n_vel + i];
      z = momentum_beta[j*n_vel + i];

      momentum_alpha[j*n_vel + i] = y * cos_PA + z * sin_PA;
      momentum_beta[j*n_vel + i] = -y * sin_PA + z * cos_PA;

      photocenter_alpha[j*n_vel + i] = momentum_alpha[j*n_vel + i] / (fline[j*n_vel + i] + EPS);
      photocenter_beta[j*n_vel + i]  = momentum_beta[j*n_vel + i] / (fline[j*n_vel + i] + EPS);

      flux_ratio = fline[j*n_vel + i]/fcon;
      ratio = flux_ratio/(1.0+flux_ratio) / DA;
      for(m=0; m<n_base; m++)
      {
        phase[j*n_vel*n_base + m*n_vel + i] =-( base[j*n_base*2 + m*2 + 0] * photocenter_alpha[j*n_vel + i]
                                               +base[j*n_base*2 + m*2 + 1] * photocenter_beta[j*n_vel + i]) * ratio;
      }
    }
  }
  
  sarm_smooth_run(pm, vel_sa, fline, n_vel, n_line, phase, n_base);
  return;
}

/*
 * smooth the profile and phase with the instrument broadening funcition.
 *
 */
int nd_fft_sarm, nd_fft_sarm_cal, npad_sarm;

fftw_complex *sarm_data_fft, *sarm_resp_fft0, *sarm_resp_fft, *sarm_conv_fft;
fftw_plan sarm_pdata, sarm_presp, sarm_pback;
double *sarm_real_data, *sarm_real_resp, *sarm_real_conv;

void sarm_smooth_init(int n_v_sarm, const double *v_sarm, double sigV)
{
  int i;
  double dV;

  npad_sarm = fmin(n_v_sarm * 0.1, 20);
  nd_fft_sarm = n_v_sarm+npad_sarm*2;
  nd_fft_sarm_cal = nd_fft_sarm/2+1;

  sarm_data_fft = (fftw_complex *) fftw_malloc((nd_fft_sarm_cal) * sizeof(fftw_complex));
  sarm_resp_fft0 = (fftw_complex *) fftw_malloc((nd_fft_sarm_cal) * sizeof(fftw_complex));
  sarm_resp_fft = (fftw_complex *) fftw_malloc((nd_fft_sarm_cal) * sizeof(fftw_complex));
  sarm_conv_fft = (fftw_complex *) fftw_malloc((nd_fft_sarm_cal) * sizeof(fftw_complex));

  sarm_real_data = (double *)fftw_malloc(nd_fft_sarm * sizeof(double));
  //real_resp = (double *)fftw_malloc(nd_fft_sa * sizeof(double));
  sarm_real_conv = (double *)fftw_malloc(nd_fft_sarm * sizeof(double));

  sarm_pdata = fftw_plan_dft_r2c_1d(nd_fft_sarm, sarm_real_data, sarm_data_fft, FFTW_PATIENT);
  //presp = fftw_plan_dft_r2c_1d(nd_fft_sa, sa_real_resp, sa_resp_fft, FFTW_PATIENT);
  sarm_pback = fftw_plan_dft_c2r_1d(nd_fft_sarm, sarm_conv_fft, sarm_real_conv, FFTW_PATIENT);

  dV = v_sarm[1] - v_sarm[0];
  for(i=0; i<nd_fft_sarm_cal; i++)
  {
    sarm_resp_fft0[i][0] = exp(-2.0 * PI*PI * sigV/dV*sigV/dV * i*i*1.0/nd_fft_sarm/nd_fft_sarm)/nd_fft_sarm;
    sarm_resp_fft0[i][1] = 0.0;
  }
  return;
}

void sarm_smooth_end()
{
  fftw_destroy_plan(sarm_pdata);
  fftw_destroy_plan(sarm_pback);

  fftw_free(sarm_data_fft);
  fftw_free(sarm_resp_fft);
  fftw_free(sarm_resp_fft0);
  fftw_free(sarm_conv_fft);
  
  fftw_free(sarm_real_data);
  //fftw_free(sarm_real_resp);
  fftw_free(sarm_real_conv);
  return;
}

void sarm_smooth_run(const void *pm, double *v_sarm, double *F_sarm, int n_v_sarm, int n_line_sarm, double *p_sarm, int n_base_sarm)
{
  int i, j, k;
  double *pmodel = (double *)pm;
  double CO, dV;

  /* 
   * line center offset: V+dV = (w-(w0+dw0))/(w0+dw0) ==> dV = - dw0/w0 
   * equivalent to redshift offset:  dz = -(1+z) dw0/w0
   */       
  CO = -pmodel[num_params_blr + num_params_sa_blr_model+3]/parset.sa_linecenter * C_Unit; 
  dV = v_sarm[1] - v_sarm[0];
  for(i=0; i<nd_fft_sarm_cal; i++)
  {
    sarm_resp_fft[i][0] = sarm_resp_fft0[i][0] * cos(2.0*PI*CO/dV * i*1.0/nd_fft_sarm);
    sarm_resp_fft[i][1] =-sarm_resp_fft0[i][0] * sin(2.0*PI*CO/dV * i*1.0/nd_fft_sarm);
  }

  /* first profile */
  for(k=0; k<n_line_sarm; k++)
  {
    memcpy(sarm_real_data+npad_sarm, F_sarm+k*n_v_sarm, n_v_sarm*sizeof(double));
    for(i=0; i<npad_sarm; i++)
      sarm_real_data[i] = sarm_real_data[nd_fft_sarm-1-i] = 0.0;
    
    /* FFT of line */
    fftw_execute(sarm_pdata);
    /* complex multiply and inverse FFT 
     * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
     * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
     */
    for(i=0; i<nd_fft_sarm_cal; i++)
    {
      sarm_conv_fft[i][0] = sarm_data_fft[i][0]*sarm_resp_fft[i][0] - sarm_data_fft[i][1]*sarm_resp_fft[i][1];
      sarm_conv_fft[i][1] = sarm_data_fft[i][0]*sarm_resp_fft[i][1] + sarm_data_fft[i][1]*sarm_resp_fft[i][0];
    }
    fftw_execute(sarm_pback);
    memcpy(F_sarm+k*n_v_sarm, sarm_real_conv+npad_sarm, n_v_sarm*sizeof(double));
  }

  /* then phase */
  for(k=0; k<n_line_sarm; k++)
  {
    for(j=0; j<n_base_sarm; j++)
    {
      memcpy(sarm_real_data+npad_sarm, p_sarm+k*n_v_sarm*n_base_sarm + j*n_v_sarm, n_v_sarm*sizeof(double));
      for(i=0; i<npad_sarm; i++)
        sarm_real_data[i] = sarm_real_data[nd_fft_sarm-1-i] = 0.0;
      
      /* FFT of line */
      fftw_execute(sarm_pdata);
      /* complex multiply and inverse FFT 
       * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
       * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
       */
      for(i=0; i<nd_fft_sarm_cal; i++)
      {
        sarm_conv_fft[i][0] = sarm_data_fft[i][0]*sarm_resp_fft[i][0] - sarm_data_fft[i][1]*sarm_resp_fft[i][1];
        sarm_conv_fft[i][1] = sarm_data_fft[i][0]*sarm_resp_fft[i][1] + sarm_data_fft[i][1]*sarm_resp_fft[i][0];
      }
      fftw_execute(sarm_pback);
      memcpy(p_sarm+k*n_v_sarm*n_base_sarm + j*n_v_sarm, sarm_real_conv+npad_sarm, n_v_sarm*sizeof(double));
    }
  }
}
#endif