/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file transfun.c
 *  \brief calculate transfer functions and emission lines
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_interp.h>

#include "brains.h"

/* 
 * interpolating continuum light curve at time tp
 */
inline double interp_con_rm(double tp)
{
  if(tp < Tcon[0])
  {
    return Fcon_rm[0];
  }
  else if(tp < Tcon[parset.n_con_recon-1])
  {
    return gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, tp, gsl_acc);
  }
  else 
  {
    return Fcon_rm[parset.n_con_recon-1];
  }
}

/* 
 * calculate A*(Fcon + ftrend)
 * 
 */
void calculate_con_rm(const void *pm)
{
  int i, m;
  double fcon, A, fcon_mean, ftrend, a0=0.0, tmp;
  double *pmodel = (double *)pm;

  A=exp(pmodel[idx_resp]);
  fcon_mean=pmodel[idx_resp + 1];

  /* add different trend in continuum and emission */
  if(parset.flag_trend_diff > 0)
  {
    tmp = 0.0;
    for(m=1; m<num_params_difftrend+1; m++)
    {
      tmp += pmodel[idx_difftrend + m-1] * pow_Tcon_data[m-1];
    }
    a0 = -tmp;

    for(i=0; i<parset.n_con_recon; i++)
    {
      ftrend = a0;
      tmp = 1.0;
      for(m=1; m<num_params_difftrend+1; m++)
      {
        tmp *= (Tcon[i] - Tmed_data);
        ftrend += pmodel[idx_difftrend + m-1] * tmp;
      }
      
      fcon = Fcon[i] + ftrend;
      /* make sure fcon is always positive */
      if(fcon < 0.0)fcon = 0.0;
      Fcon_rm[i] = A * (fcon - fcon_mean);
    }
  }
  else 
  {
    for(i=0; i<parset.n_con_recon; i++)
    {
      /* if fcon is negative, its pow raises an error */
      Fcon_rm[i] = A * (Fcon[i] - fcon_mean);
    }
  }

  return;
}

/*!
 * This function calculate 1d line from a given BLR model.
 *
 * Note that continuum light curve has been obtained in advance 
 *
 * Input:
 *   pm: model parameters
 *   Tl: epochs of line
 *   nl: number of points
 * Output:
 *   Fl: line fluxes
 *
 */
void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl)
{
  int i, j;
  double A, fcon_mean, fline, fline_mean, fcon_rm, tl, tc, tau, dTransTau;
  double *pmodel = (double *)pm;
  A=exp(pmodel[idx_resp]);
  fcon_mean=pmodel[idx_resp + 1];
  
  fline_mean = Fline_mean * A * fcon_mean; /* scale with the response coefficient */

  dTransTau = TransTau[1] - TransTau[0];

  for(i=0;i<nl;i++)
  {
    tl = Tl[i];
    fline = 0.0;
    for(j=0; j<parset.n_tau; j++)
    {
      tau = TransTau[j];
      tc = tl - tau;
      //fcon_rm = gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, tc, gsl_acc); /* interpolation */
      fcon_rm = interp_con_rm(tc);
    
      fline += Trans1D[j] * fcon_rm;     /*  line response */
    }
    Fl[i] = fline * dTransTau + fline_mean;  /* add back mean spectrum */
  }

  return;
}

/*!
 * This function caclulate 2d line from obtained transfer function.
 */
void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv)
{
  int i, j, k;
  double A, fcon_mean, tau, tl, tc, fcon_rm, fnarrow, dTransTau, sigV_mean, linecenter_mean;
  double *pmodel = (double *)pm;
  A=exp(pmodel[idx_resp]);
  fcon_mean=pmodel[idx_resp + 1];

  /* multiply response coefficient */
  for(i=0; i<nv; i++)
  {
    Fline2d_mean_buf[i] = Fline2d_mean[i] * A * fcon_mean; /* scale with the response coefficient */
  }

  dTransTau = TransTau[1] - TransTau[0];

  for(j=0;j<nl; j++)
  {
    for(i=0; i<nv; i++)
      fl2d[j*nv + i] = 0.0;

    tl = Tl[j];
    for(k=0; k<parset.n_tau; k++)
    {
      tau = TransTau[k];
      tc = tl - tau;
      //fcon_rm = gsl_interp_eval(gsl_linear, Tcon, Fcon_rm, tc, gsl_acc);
      fcon_rm = interp_con_rm(tc);

      for(i=0; i<nv; i++)
      {
        fl2d[j*nv + i] += trans2d[k*nv+i] * fcon_rm;
      }
    }
    for(i=0; i<nv; i++)
    {
      fl2d[j*nv + i] *= dTransTau;
      fl2d[j*nv + i] += Fline2d_mean_buf[i]; /* add back mean spectrum */
    }
  }

  /* add intrinsic narrow line */
  if(parset.flag_narrowline != 0)
  {
    double flux, width, shift;
    if(parset.flag_narrowline == 1)  /* fixed narrow line */
    {
      flux = parset.flux_narrowline;
      width = parset.width_narrowline;
      shift = parset.shift_narrowline;
    }
    else if(parset.flag_narrowline == 2) /* narrow line with Gaussian priors */
    {
      flux =  parset.flux_narrowline  + pmodel[num_params_blr_model] * parset.flux_narrowline_err;
      width = parset.width_narrowline + pmodel[num_params_blr_model+1] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr_model+2] * parset.shift_narrowline_err;
    }
    else  /* narrow line with logrithmic prior of flux */
    {
      flux =  exp(pmodel[num_params_blr_model]);
      width = parset.width_narrowline + pmodel[num_params_blr_model+1] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr_model+2] * parset.shift_narrowline_err;
    }

    width = fmax(1.0e-10, width); /* make sure thant width is not zero */

    for(i=0; i<nv; i++)
    {
      fnarrow = flux/sqrt(2.0*PI)/width * exp( -0.5 * pow( (transv[i] - shift)/(width), 2.0) );
      for(j = 0; j<nl; j++)
      {
        fl2d[j*nv + i] += fnarrow;
      }
    } 
  }

  /* smooth the line profile, get the mean sigV and linecenter */
  line_gaussian_smooth_2D_FFT(transv, fl2d, nl, nv, pm);
  /* smooth the mean line profile */
  // line_gaussian_smooth_FFT_width_linecenter(transv, Fline2d_mean_buf, nv, sigV_mean, linecenter_mean);
}

/*!
 * This function calculates 1d transfer function.
 */
void transfun_1d_cal_cloud(const void *pm, int flag_save)
{
  /* generate cloud sample and calculate the corresponding time lags and weights */
  gen_cloud_sample(pm, 1, flag_save);
  transfun_1d_cal_with_sample();
  return;
}

/*!
 * This function calculates 1d transfer function.
 */
void transfun_1d_cal_with_sample()
{
  int i, idt;
  double tau_min, tau_max, dTransTau;
  double Anorm, dis, fline_mean;

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
    Trans1D[i] = 0.0;
  }
  
  fline_mean = 0.0;
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    dis = clouds_tau[i];
    idt = (dis - tau_min)/dTransTau;
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += clouds_weight[i];

    fline_mean += clouds_weight_mean[i];
  }

  /* normalize transfer function */
  Anorm = fline_mean * dTransTau + EPS;

  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] /= Anorm;
  }
  
  Fline_mean = 1.0; 

  smooth_transfer_function_tau(Trans1D);
  
  return;
}


/*!
 * This function calculates 2d transfer function.
 */
void transfun_2d_cal_cloud(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  /* generate cloud sample and calculate the corresponding time lags, LOS velocity, and weights */
  gen_cloud_sample(pm, 2, flag_save);
  transfun_2d_cal_with_sample(transv, trans2d, n_vel);
  return;
}

/*!
 * This function calculates 2d transfer function.
 */
void transfun_2d_cal_with_sample(double *transv, double *trans2d, int n_vel)
{
  int i, j, idt, idV;
  double tau_min, tau_max, dTransTau;
  double Anorm, dis, V, V_offset, dV;

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
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   /* cleanup of transfer function */
  
  for(j=0; j<n_vel; j++)
  {
    Fline2d_mean[j] = 0.0;
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

      Fline2d_mean[idV] += clouds_weight_mean[i];
    }
  }

  /* normalize transfer function */
  Anorm = EPS;
  for(j=0; j<n_vel; j++)
  {
    Anorm += Fline2d_mean[j];
  }
  Anorm *= dV;

  for(j=0; j<n_vel; j++)
  {
    Fline2d_mean[j] /= Anorm;
  }
  
  Anorm *=  dTransTau;
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<n_vel; j++)
    {
      trans2d[i*n_vel+j] /= Anorm;
    }
  }

  smooth_transfer_function2d_tau(trans2d, n_vel);

  return;
}

/* 
 * smooth the transfer function along time lag axis
 * this is to surppress statistic noises
 */
void smooth_transfer_function_tau(double *trans1d)
{
  memcpy(hist_in->data, trans1d, parset.n_tau*sizeof(double));
  gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 0, hist_in, hist_out, gauss_p);
  memcpy(trans1d, hist_out->data, parset.n_tau*sizeof(double));

  return;
}
void smooth_transfer_function2d_tau(double *trans2d, int n_vel)
{
  int i, j;
  for(j=0; j<n_vel; j++)
  {
    for(i=0; i<parset.n_tau; i++)
    {
      hist_in->data[i] = trans2d[i*n_vel + j];
    }
    gsl_filter_gaussian(GSL_FILTER_END_PADVALUE, alpha, 0, hist_in, hist_out, gauss_p);
    for(i=0; i<parset.n_tau; i++)
    {
      trans2d[i*n_vel + j] = hist_out->data[i];
    }
  }
  return;
}

/*!
 * This function calculates emissivity-weighted 2d transfer function.
 */
void transfun_2d_ew_cal_cloud(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  /* generate cloud sample and calculate the corresponding time lags, LOS velocity, and weights */
  gen_cloud_sample(pm, 2, flag_save);
  transfun_2d_ew_cal_with_sample(transv, trans2d, n_vel);
  return;
}
/*!
 * This function calculates emissivity-weighted transfer function.
 *
 */
void transfun_2d_ew_cal_with_sample(double *transv, double *trans2d, int n_vel)
{
  int i, j, idt, idV;
  double tau_min, tau_max, dTransTau;
  double Anorm, dis, V, V_offset, dV;

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
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   /* cleanup of transfer function */

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
      trans2d[idt*n_vel + idV] += clouds_weight_mean[i];
    }
  }

  /* normalize transfer function */
  Anorm = EPS;
  for(i=0; i<parset.n_tau; i++)
    for(j=0; j<n_vel; j++)
      Anorm += trans2d[i*n_vel + j];
    
  Anorm *=  dTransTau * dV;
  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<n_vel; j++)
    {
      trans2d[i*n_vel+j] /= Anorm;
    }
  }

  smooth_transfer_function2d_tau(trans2d, n_vel);

  return;
}

/*
 * calculate centroid time lag of 1d transfer function
 *
 */
double calculate_cent_transfun_1d(double *transtau, double *trans1d, int ntau)
{
  int j;
  double norm, r_tau;
  // calculate centroid time lag
  norm = 0.0;
  r_tau = 0.0;
  for(j=0; j<ntau; j++)
  {
    norm += trans1d[j];
    r_tau += trans1d[j] * transtau[j];
  }
  r_tau = r_tau/norm;

  return r_tau;
}