/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
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
  int i, j, m;
  double fline, fcon, tl, tc, tau, A, Ag, ftrend, a0=0.0, tmp;
  double *pmodel = (double *)pm;

  A=exp(pmodel[num_params_blr-3]); /*  response coefficient */
  Ag=pmodel[num_params_blr-2];     /*  no-linearity of response */

  if(parset.flag_trend_diff > 0)
  {
    tmp = 0.0;
    for(m=1; m<num_params_difftrend+1; m++)
    {
      tmp += pmodel[num_params_blr + num_params_drw + num_params_trend + m-1] * pow_Tcon_data[m-1];
    }
    a0 = -tmp;
  }

  for(i=0;i<nl;i++)
  {
    tl = Tl[i];
    fline = 0.0;
    for(j=0; j<parset.n_tau; j++)
    {
      tau = TransTau[j];
      tc = tl - tau;
      if(tc>=Tcon_min && tc <=Tcon_max)
      {
        fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, tc, gsl_acc); /* interpolation */
      }
      else
      {
        /* fcon = mean */ 
        fcon = con_q[0];
        tmp = 1.0;
        for(m=1; m < nq; m++)/*  beyond the range, set to be the long-term trend */
        {
          tmp *= tc;
          fcon += con_q[m] * tmp;
        }
      }
      
      /* add different trend in continuum and emission */
      if(parset.flag_trend_diff > 0)
      {
        ftrend = a0;
        tmp = 1.0;
        for(m=1; m<num_params_difftrend+1; m++)
        {
          tmp *= (tc - Tmed_data);
          ftrend += pmodel[num_params_blr + num_params_drw + num_params_trend + m-1] * tmp;
        }
        fcon += ftrend;
      }

      /* take into account the possibility that fcon is negative */
      if(fcon > 0.0)
      {
        fline += Trans1D[j] * fcon * pow(fcon, Ag);     /*  line response */
      }	
    }
    Fl[i] = fline * dTransTau * A;
  }

  return;
}

/*!
 * This function caclulate 2d line from obtained transfer function.
 */
void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv)
{
  int i, j, k, m;
  double tau, tl, tc, fcon, A, Ag, ftrend, fnarrow, a0=0.0, tmp;
  double *pmodel = (double *)pm;

  A=exp(pmodel[num_params_blr-3]);
  Ag=pmodel[num_params_blr-2];

  if(parset.flag_trend_diff > 0)
  {
    tmp = 0.0;
    for(m=1; m<num_params_difftrend+1; m++)
    {
      tmp += pmodel[num_params_blr + num_params_drw + num_params_trend + m-1] * pow_Tcon_data[m-1];
    }
    a0 = -tmp;
  }

  for(j=0;j<nl; j++)
  {
    for(i=0; i<nv; i++)
      fl2d[j*nv + i] = 0.0;

    tl = Tl[j];
    for(k=0; k<parset.n_tau; k++)
    {
      tau = TransTau[k];
      tc = tl - tau;
      if(tc>=Tcon_min && tc <=Tcon_max)
      {
        fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, tc, gsl_acc);          
      }
      else
      {
        /* fcon = mean */ 
        fcon = con_q[0];
        tmp = 1.0;
        for(m=1; m < nq; m++)/*  beyond the range, set to be the long-term trend */
        {
          tmp *= tc;
          fcon += con_q[m] * tmp;
        }
      }

      /* add different trend in continuum and emission line */
      if(parset.flag_trend_diff > 0)
      {
        ftrend = a0;
        tmp = 1.0;
        for(m=1; m<num_params_difftrend+1; m++)
        {
          tmp *= (tc - Tmed_data);
          ftrend += pmodel[num_params_blr + num_params_drw + num_params_trend + m-1] * tmp;
        }
        fcon += ftrend;
      }

      for(i=0; i<nv; i++)
      {
        if(fcon > 0.0)
        {
           fl2d[j*nv + i] += trans2d[k*nv+i] * fcon * pow(fcon, Ag);
           //fline += trans2d[k*nv+i] * fcon;
        }
      }
    }
    for(i=0; i<nv; i++)
    {
      fl2d[j*nv + i] *= dTransTau * A;
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
      fnarrow = flux * exp( -0.5 * pow( (transv[i] - shift)/(width), 2.0) );
      for(j = 0; j<nl; j++)
      {
        fl2d[j*nv + i] += fnarrow;
      }
    } 
  }

  /* smooth the line profile */
  line_gaussian_smooth_2D_FFT(transv, fl2d, nl, nv, pm);
}

/*!
 * This function calculates 1d transfer function.
 */
void transfun_1d_cal(const void *pm, int flag_save)
{
  int i, idt;
  double Anorm, dis;

  /* generate cloud sample and calculate the corresponding time lags and weights */
  transfun_1d_cloud_sample(pm, flag_save);

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

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    dis = clouds_tau[i];
    idt = (dis - tau_min)/dTransTau;
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += clouds_weight[i];
  }

  /* normalize transfer function */
  Anorm = 0.0;
  for(i=0;i<parset.n_tau;i++)
  {
    Anorm += Trans1D[i];
  }
  Anorm *= dTransTau;
  /* check if we get a zero transfer function */
  if(Anorm > 0.0)
  {
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] /= Anorm;
    }
  }
  else
  {
    printf(" Warning, zero 1d transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }

  return;
}


/*!
 * This function calculates 2d transfer function.
 */
void transfun_2d_cal(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idt, idV;
  double Anorm, dis, V, dV;

  /* generate cloud sample and calculate the corresponding time lags, LOS velocity, and weights */
  transfun_2d_cloud_sample(pm, transv, trans2d, n_vel, flag_save);

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
      if(V<transv[0] || V >= transv[n_vel-1]+dV)
        continue;
      idV = (V - transv[0])/dV;
      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += clouds_weight[i];
    }
  }

  /* normalize transfer function */
  Anorm = 0.0;
  for(i=0; i<parset.n_tau; i++)
    for(j=0; j<n_vel; j++)
    {
      Anorm += trans2d[i*n_vel+j];
    }
  Anorm *= (dV * dTransTau);

  /* check if we get a zero transfer function */
  if(Anorm > 0.0)
  {
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] /= Anorm;
      }
    }
  }
  else
  {
    printf(" Warning, zero 2d transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }
  return;
}

