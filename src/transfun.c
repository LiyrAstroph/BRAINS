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

#include "allvars.h"
#include "proto.h"

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
  int i, j, k;
  double fline, fcon, tl, tc, tau, A, Ag, ftrend, a0, tmp;
  double *pmodel = (double *)pm;

  A=exp(pmodel[0]); /*  response coefficient */
  Ag=pmodel[1];     /*  no-linearity of response */

  if(parset.flag_trend_diff > 0)
  {
    tmp = 0.0;
    for(k=1; k<num_params_difftrend+1; k++)
    {
      tmp += pmodel[num_params_blr + 4 + parset.flag_trend + k-1] * pow_Tcon_data[k-1];
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
        for(k=1; k < nq; k++)/*  beyond the range, set to be the long-term trend */
        {
          tmp *= tc;
          fcon += con_q[k] * tmp;
        }
      }
      
      /* add different trend in continuum and emission */
      if(parset.flag_trend_diff > 0)
      {
        ftrend = a0;
        tmp = 1.0;
        for(k=1; k<num_params_difftrend+1; k++)
        {
          tmp *= (tc - Tmed_data);
          ftrend += pmodel[num_params_blr + 4 + parset.flag_trend + k-1] * tmp;
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
  double tau, tl, tc, fcon, A, Ag, ftrend, fnarrow, a0, tmp;
  double *pmodel = (double *)pm;

  A=exp(pmodel[0]);
  Ag=pmodel[1];

  if(parset.flag_trend_diff > 0)
  {
    tmp = 0.0;
    for(k=1; k<num_params_difftrend+1; k++)
    {
      tmp += pmodel[num_params_blr + 4 + parset.flag_trend + k-1] * pow_Tcon_data[k-1];
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
        for(k=1; k<num_params_difftrend+1; k++)
        {
          tmp *= (tc - Tmed_data);
          ftrend += pmodel[num_params_blr + 4 + parset.flag_trend + k-1] * tmp;
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
      flux =  parset.flux_narrowline  + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-3] * parset.flux_narrowline_err;
      width = parset.width_narrowline + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-2] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-1] * parset.shift_narrowline_err;
    }
    else  /* narrow line with logrithmic prior of flux */
    {
      flux =  exp(pmodel[num_params_blr-num_params_res-num_params_linecenter-1-3]);
      width = parset.width_narrowline + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-2] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-1] * parset.shift_narrowline_err;
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

  tau_min = tmp_tau[0];
  tau_max = tmp_tau[0];
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    if(tau_min > tmp_tau[i])
      tau_min = tmp_tau[i];
    if(tau_max < tmp_tau[i])
      tau_max = tmp_tau[i];
  }

  dTransTau = (tau_max - tau_min)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = tau_min + dTransTau * i;
    Trans1D[i] = 0.0;
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    dis = tmp_tau[i];
    idt = (dis - tau_min)/dTransTau;
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += tmp_weight[i];
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

  tau_min = tmp_tau[0];
  tau_max = tmp_tau[0];
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    if(tau_min > tmp_tau[i])
      tau_min = tmp_tau[i];
    if(tau_max < tmp_tau[i])
      tau_max = tmp_tau[i];
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
    dis = tmp_tau[i];
    idt = (dis - tau_min)/dTransTau;

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      V = tmp_vel[i*parset.n_vel_per_cloud + j];
      if(V<transv[0] || V >= transv[n_vel-1]+dV)
        continue;
      idV = (V - transv[0])/dV;
      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += tmp_weight[i];
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
}

/*================================================================
 * model 1
 * Brewer et al. (2011)'s model
 *
 * geometry: radial Gamma distribution
 * dynamics: elliptical orbits
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model1(const void *pm, int flag_save)
{
  int i, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe;
  double weight, rnd;
  double *pr;
  BLRmodel1 *model = (BLRmodel1 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  beta = model->beta;         
  F = model->F;
  mu = exp(model->mu);                 /* mean radius */
  k = model->k;  

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  /* "which_parameter_update = -1" means that all parameters are updated, 
   * usually occurs at the initial step.
   */
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

   /* right-handed framework
    * first rotate around y axis by an angle of Lthe, then roate around z axis 
    * by an angle of Lphi
    */
  
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

    /* conter-rotate around y, LOS is x-axis */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }
  
  return;
}


/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_sample_model1(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, L, E, linecenter = 0.0;
  double V, weight, rnd;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel1 *model = (BLRmodel1 *)pm;
  double Emin, Lmax, Vr, Vr2, Vph, mbh, chi, lambda, q;
  

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = acos(model->inc);
  beta = model->beta;
  F = model->F;
  mu = exp(model->mu);
  k = model->k;

  a = 1.0/beta/beta;
  s = mu/a;
  rin = mu*F;
  sig = (1.0-F)*s;

  mbh = exp(model->mbh);
  lambda = model->lambda;
  q = model->q;
  
  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

  /*xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =  sin(Lthe) * x + cos(Lthe) * z;*/

    xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb =  sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

    /* conter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;
    
    /* velocity  
     * note that a cloud moves in its orbit plane, whose direction
     * is determined by the direction of its angular momentum.
     */
    Emin = - mbh/r;
    //Ecirc = 0.5 * Emin;
    //Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      chi = lambda * gsl_ran_ugaussian(gsl_r);
      E = Emin / (1.0 + exp(-chi));

      Lmax = sqrt(2.0 * r*r * (E + mbh/r));      

      if(lambda>1.0e-2)   /* make sure that exp is caculatable. */
        L = Lmax * lambda * log( (exp(1.0/lambda) - 1.0) * gsl_rng_uniform(gsl_r) + 1.0 );
      else
        L = Lmax * (1.0 + lambda * log(gsl_rng_uniform(gsl_r)) );
 
      Vr2 = 2.0 * (E + mbh/r) - L*L/r/r;
      if(Vr2>=0.0)
      {
        u = gsl_rng_uniform(gsl_r);
        Vr = sqrt(Vr2) * (u<q?-1.0:1.0);
      }
      else
      {
        Vr = 0.0;
      }
      Vph = L/r; /* RM cannot distiguish the orientation of the rotation. */

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if(zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      
      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}



/*================================================================
 * model 2
 * 
 * geometry: radial Gamma distribution
 * dynamics: ellipitical orbits (Gaussian around circular orbits)
 *================================================================
 */
/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_sample_model2(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, linecenter = 0.0;
  double V, weight, rnd;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel2 *model = (BLRmodel2 *)pm;
  double Emin, Ecirc, Lcirc, Vcirc, Vr, Vph, mbh, sigr, sigtheta, rhor, rhotheta;
  

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = acos(model->inc);
  beta = model->beta;
  F = model->F;
  mu = exp(model->mu);
  k = model->k;

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  mbh = exp(model->mbh);
  sigr = model->sigr;
  sigtheta = model->sigtheta * PI;
  
  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }


  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum */    
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;


  /*xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =  sin(Lthe) * x + cos(Lthe) * z;*/

    xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb =  sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

    /* conter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;
    
    
    /* velocity  
     * note that a cloud moves in its orbit plane, whose direction
     * is determined by the direction of its angular momentum.
     */
    Emin = - mbh/r;
    Ecirc = 0.5 * Emin;
    Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
    Vcirc = Lcirc/r;
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {    
      rhor = gsl_ran_ugaussian(gsl_r) * sigr + 1.0;
      rhotheta = gsl_ran_ugaussian(gsl_r) * sigtheta + 0.5*PI;

      Vr = rhor * cos(rhotheta) * Vcirc;
      Vph = rhor * sin(rhotheta) * Vcirc; //RM cannot distiguish the orientation of the rotation.

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if(zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);


      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}


/*====================================================================
 * model 3
 *
 * geometry: radial power law
 * dynamics: 
 *====================================================================
 */

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model3(const void *pm, int flag_save)
{
  int i, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe;
  double weight, rnd;
  double *pr;
  BLRmodel3 *model = (BLRmodel3 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  alpha = model->alpha;         
  F = exp(model->F);
  Rin = exp(model->Rin);                 /* mean radius */
  k = model->k;  
  
  if(F*Rin > rcloud_max_set)
    F = rcloud_max_set/Rin;

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  
  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));
    
    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

    /* right-handed framework
     * first rotate around y axis by an angle of Lthe, then roate around z axis 
     * by an angle of Lphi
     */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z;*/

    xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

    /* conter-rotate around y, LOS is x-axis */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  return;
}

/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_sample_model3(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, linecenter=0.0;
  double V, weight, rnd;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel3 *model = (BLRmodel3 *)pm;
  double Emin, Lmax, Vr, Vph, mbh, xi, q;
  

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = acos(model->inc);
  alpha = model->alpha;
  F = model->F;
  Rin = exp(model->Rin);
  k = model->k;

  mbh = exp(model->mbh);
  xi = model->xi;
  q = model->q;
  
  if(F*Rin > rcloud_max_set)
    F = rcloud_max_set/Rin;

  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;


  /*xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =  sin(Lthe) * x + cos(Lthe) * z;*/

    xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb =  sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

// conter-rotate around y
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;
    
// velocity  
// note that a cloud moves in its orbit plane, whose direction
// is determined by the direction of its angular momentum.
    Emin = - mbh/r;
    //Ecirc = 0.5 * Emin;
    //Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      
      E = Emin / 2.0;
      Lmax = sqrt(2.0 * r*r * (E + mbh/r));      
      L = Lmax;
 
      Vr = xi * sqrt(2.0*mbh/r);
      if(q <= 0.5)
        Vr = -Vr;
      Vph = L/r; //RM cannot distiguish the orientation of the rotation.

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if(zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}

/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_sample_model4(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, linecenter=0.0;
  double V, weight, rnd;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel4 *model = (BLRmodel4 *)pm;
  double Emin, Lmax, Vr, Vph, mbh, xi, q;
  

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = acos(model->inc);
  alpha = model->alpha;
  F = model->F;
  Rin = exp(model->Rin);
  k = model->k;

  mbh = exp(model->mbh);
  xi = model->xi;
  q = model->q;

  if(Rin*F > rcloud_max_set)
    F = rcloud_max_set/Rin;
  
  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }
  
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum   */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;


  /*xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =  sin(Lthe) * x + cos(Lthe) * z; */
    
    xb =  cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb =  sin(Lthe) * x;

    zb0 = zb;
    if(zb0 < 0.0)
      zb = -zb;

    /* conter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;
    
// velocity  
// note that a cloud moves in its orbit plane, whose direction
// is determined by the direction of its angular momentum.
    Emin = - mbh/r;
    //Ecirc = 0.5 * Emin;
    //Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      
      E = Emin / 2.0;
      Lmax = sqrt(2.0 * r*r * (E + mbh/r));      
      L = Lmax;
 
      Vr = xi * sqrt(2.0*mbh/r);
      if(q<=0.5)
        Vr = -Vr;
      Vph = sqrt(1.0-2.0*xi*xi) * L/r; //RM cannot distiguish the orientation of the rotation.

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if(zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}

/*====================================================================
 * model 5
 *
 * geometry: double pow-law 
 * dynamics: ellipitical orbits and inflow/outflow as in Pancoast's model
 *====================================================================
 */

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model5(const void *pm, int flag_save)
{
  int i, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double Lphi, Lthe, cos_Lphi, sin_Lphi, cos_Lthe, sin_Lthe, cos_inc_cmp, sin_inc_cmp;
  double weight, rndr, rnd, rnd_xi, rnd_frac, frac1, frac2, ratio;
  double *pr;
  BLRmodel5 *model = (BLRmodel5 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  alpha = model->alpha;         
  Fin = model->Fin;
  Fout = exp(model->Fout);                 
  mu = exp(model->mu);            /* mean radius */
  k = model->k;  
  gam = model-> gam;
  xi = model->xi;
  

  if(mu*Fout > rcloud_max_set) // cloud radius samller than rcloud_max_set
    Fout = rcloud_max_set/mu;

  frac1 = 1.0/(alpha+1.0) * (1.0 - pow(Fin, alpha+1.0));
  frac2 = 1.0/(alpha-1.0) * (1.0 - pow(Fout, -alpha+1.0));
  ratio = frac1/(frac1 + frac2);

  sin_inc_cmp = cos(inc);//sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc);//cos(PI/2.0 - inc);

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        
        rnd_frac = gsl_rng_uniform(gsl_r);
        rnd = gsl_rng_uniform(gsl_r);
        if(rnd_frac < ratio)
        {
          rndr = pow( 1.0 - rnd * (1.0 - pow(Fin, alpha+1.0)), 1.0/(1.0+alpha));
        }
        else
        {
          rndr = pow( 1.0 - rnd * (1.0 - pow(Fout, -alpha+1.0)), 1.0/(1.0-alpha));
        }
        r = rndr*mu;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    cos_Lphi = cos(Lphi);
    sin_Lphi = sin(Lphi);
    cos_Lthe = cos(Lthe);
    sin_Lthe = sin(Lthe);
    
    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z;*/

    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;
    
    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  return;
}

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_2d_cloud_sample_model5(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, V, Vr, Vph, Vkep, rhoV, theV, linecenter=0.0;
  double cos_Lphi, sin_Lphi, cos_Lthe, sin_Lthe, cos_inc_cmp, sin_inc_cmp;
  double weight, rndr, rnd, rnd_frac, rnd_xi, frac1, frac2, ratio, Rs, g;
  double vx, vy, vz, vxb, vyb, vzb;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel5 *model = (BLRmodel5 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  alpha = model->alpha;         
  Fin = model->Fin;
  Fout = exp(model->Fout);                 
  mu = exp(model->mu);            /* mean radius */
  k = model->k;  
  gam = model->gam;
  xi = model->xi;

  mbh = exp(model->mbh);
  fellip = model->fellip;
  fflow = model->fflow;
  sigr_circ = exp(model->sigr_circ);
  sigthe_circ = exp(model->sigthe_circ);
  sigr_rad = exp(model->sigr_rad);
  sigthe_rad = exp(model->sigthe_rad);
  theta_rot = model->theta_rot*PI/180.0;
  sig_turb = exp(model->sig_turb);

  Rs = 3.0e11*mbh/CM_PER_LD; // Schwarzchild radius in a unit of light-days

  if(mu*Fout > rcloud_max_set)
    Fout = rcloud_max_set/mu;
  
  frac1 = 1.0/(alpha+1.0) * (1.0 - pow(Fin,   alpha+1.0));
  frac2 = 1.0/(alpha-1.0) * (1.0 - pow(Fout, -alpha+1.0));
  ratio = frac1/(frac1 + frac2);

  sin_inc_cmp = cos(inc);//sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc);//cos(PI/2.0 - inc);

  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  
  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        
        rnd_frac = gsl_rng_uniform(gsl_r);
        rnd = gsl_rng_uniform(gsl_r);
        if(rnd_frac < ratio)
        {
          rndr = pow( 1.0 - rnd * (1.0 - pow(Fin, alpha+1.0)), 1.0/(1.0+alpha));
        }
        else
        {
          rndr = pow( 1.0 - rnd * (1.0 - pow(Fout, -alpha+1.0)), 1.0/(1.0-alpha));
        }
        r = rndr*mu;
        nc++;
      }

      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    cos_Lphi = cos(Lphi);
    sin_Lphi = sin(Lphi);
    cos_Lthe = cos(Lthe);
    sin_Lthe = sin(Lthe);
    
    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z;*/

    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;
    
    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// counter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    Vkep = sqrt(mbh/r);

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV = (gsl_ran_ugaussian(gsl_r) * sigthe_circ + 0.5) * PI;
      }
      else
      {
        if(fflow <= 0.5)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad + 1.0) * PI + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad ) * PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos_phi - Vph * sin_phi;
      vy = Vr * sin_phi + Vph * cos_phi;
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos_Lthe*cos_Lphi * vx + sin_Lphi * vy;
      vyb =-cos_Lthe*sin_Lphi * vx + cos_Lphi * vy;
      vzb = sin_Lthe * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
      vy = vyb;
      vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += gsl_ran_ugaussian(gsl_r) * sig_turb * Vkep; // add turbulence velocity

      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = sqrt( (1.0 + V/C_Unit) / (1.0 - V/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}


/*================================================================
 * model 6
 * Pancoast et al. (2014)'s model
 *
 * geometry: radial Gamma distribution
 * dynamics: ellipitcal orbits and inflow/outflow 
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model6(const void *pm, int flag_save)
{
  int i, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp;
  double weight, rnd, rnd_xi;
  double *pr;
  BLRmodel6 *model = (BLRmodel6 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  beta = model->beta;         
  F = model->F;
  mu = exp(model->mu);                 /* mean radius */
  k = model->k;  
  gam = model-> gam;
  xi = model->xi;

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  sin_inc_cmp = cos(inc);//sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc);//cos(PI/2.0 - inc);

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    /*x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc); */

    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);    
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  return;
}

/* 
 * This function caclulate 2d transfer function.
 */
void transfun_2d_cloud_sample_model6(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double V, rhoV, theV, Vr, Vph, Vkep, Rs, g;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp, linecenter=0.0;
  double weight, rnd, rnd_xi;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel6 *model = (BLRmodel6 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  beta = model->beta;         
  F = model->F;
  mu = exp(model->mu);                 /* mean radius */
  k = model->k;  
  gam = model-> gam;
  xi = model->xi;

  mbh = exp(model->mbh);
  fellip = model->fellip;
  fflow = model->fflow;
  sigr_circ = exp(model->sigr_circ);
  sigthe_circ = exp(model->sigthe_circ);
  sigr_rad = exp(model->sigr_rad);
  sigthe_rad = exp(model->sigthe_rad);
  theta_rot = model->theta_rot*PI/180.0;
  sig_turb = exp(model->sig_turb);

  Rs = 3.0e11*mbh/CM_PER_LD; // Schwarzchild radius in a unit of light-days

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F + Rs;  // include Scharzschild radius
  sig=(1.0-F)*s;
  
  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);

  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the 
  // initial step.
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  
  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    /* x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
       y = yb;
       z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc); */
    
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    Vkep = sqrt(mbh/r);
    
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV =  (gsl_ran_ugaussian(gsl_r) * sigthe_circ + 0.5)*PI;
      }
      else
      {
        if(fflow <= 0.5) /* inflow */
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad   + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad + 1.0) * PI + theta_rot;
        }
        else     /* outflow */
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad) * PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos_phi - Vph * sin_phi;
      vy = Vr * sin_phi + Vph * cos_phi;
      vz = 0.0;    
      
    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos_Lthe*cos_Lphi * vx + sin_Lphi * vy;
      vyb =-cos_Lthe*sin_Lphi * vx + cos_Lphi * vy;
      vzb = sin_Lthe * vx;


      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      /*vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);*/

      vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
      vy = vyb;
      vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += gsl_ran_ugaussian(gsl_r) * sig_turb * Vkep; // add turbulence velocity

      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = sqrt( (1.0 + V/C_Unit) / (1.0 - V/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}

/*================================================================
 * model 7
 * shadowed model
 *
 * geometry: radial Gamma distribution
 * dynamics: elliptical orbits and inflow/outflow as in Pancoast'model
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model7(const void *pm, int flag_save)
{
  int i, nc, flag_update=0, num_sh;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp;
  double weight, rnd, rnd_xi;
  double *pr;
  BLRmodel7 *model = (BLRmodel7 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  beta = model->beta;         
  F = model->F;
  mu = exp(model->mu);                 /* mean radius */
  k = model->k;  
  gam = model-> gam;
  xi = model->xi;

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the 
  // initial step.
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  
  /* number of particles in first region */
  
  num_sh = (int)(parset.n_cloud_per_task * model->fsh);

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<num_sh; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<num_sh; i++)
  {
    /* generate a direction of the angular momentum of the orbit */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  /* second BLR region */
  rin = exp(model->mu_un) * model->F_un;
  a = 1.0/(model->beta_un * model->beta_un);
  sig = exp(model->mu_un) * (1.0-model->F_un)/a;
  double Lopn_cos_un1, Lopn_cos_un2;
  Lopn_cos_un1 = Lopn_cos;
  if(model->opn + model->opn_un < 90.0)
    Lopn_cos_un2 = cos((model->opn + model->opn_un)/180.0*PI);
  else
    Lopn_cos_un2 = 0.0;

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i = num_sh; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=num_sh; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos_un2 + (Lopn_cos_un1-Lopn_cos_un2) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  return;
}

/* 
 * This function caclulate 2d transfer function.
 */
void transfun_2d_cloud_sample_model7(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc, flag_update=0, num_sh;
  double r, phi, dis, Lopn_cos, cos_phi, sin_phi;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb, Rs, g, sig_turb;
  double V, rhoV, theV, Vr, Vph, Vkep;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp,linecenter=0.0;
  double weight, rnd, rnd_xi;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel7 *model = (BLRmodel7 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  beta = model->beta;         
  F = model->F;
  mu = exp(model->mu);                 /* mean radius */
  k = model->k;  
  gam = model-> gam;
  xi = model->xi;

  mbh = exp(model->mbh);
  fellip = model->fellip;
  fflow = model->fflow;
  sigr_circ = exp(model->sigr_circ);
  sigthe_circ = exp(model->sigthe_circ);
  sigr_rad = exp(model->sigr_rad);
  sigthe_rad = exp(model->sigthe_rad);
  theta_rot = model->theta_rot*PI/180.0;
  sig_turb = exp(model->sig_turb);

  Rs = 3.0e11*mbh/CM_PER_LD; // Schwarzchild radius in a unit of light-days

  a = 1.0/beta/beta;
  s = mu/a;
  rin=Rs + mu*F;
  sig=(1.0-F)*s;
  
  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);

  if(parset.flag_linecenter !=0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err; 
  }

  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  else
  {
    for(i=0; i<num_params_radial_samp;i++)
    {
      if(which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }
  
  num_sh = (int)(parset.n_cloud_per_task * model->fsh);

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=0; i<num_sh; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=0; i<num_sh; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    Vkep = sqrt(mbh/r);

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV = (gsl_ran_ugaussian(gsl_r) * sigthe_circ + 0.5) * PI;
      }
      else
      {
        if(fflow <= 0.5)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad + 1.0) *PI + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad) * PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos_phi - Vph * sin_phi;
      vy = Vr * sin_phi + Vph * cos_phi;
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos_Lthe*cos_Lphi * vx + sin_Lphi * vy;
      vyb =-cos_Lthe*sin_Lphi * vx + cos_Lphi * vy;
      vzb = sin_Lthe * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
      vy = vyb;
      vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += gsl_ran_ugaussian(gsl_r) * sig_turb * Vkep;

      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = sqrt( (1.0 + V/C_Unit) / (1.0 - V/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;

      V += linecenter;
      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  //second BLR region
  fellip = model->fellip_un;
  fflow = model->fflow_un;
  rin = exp(model->mu_un) * model->F_un + Rs;
  a = 1.0/(model->beta_un * model->beta_un);
  sig = exp(model->mu_un) * (1.0-model->F_un)/a;
  double Lopn_cos_un1, Lopn_cos_un2;
  Lopn_cos_un1 = Lopn_cos;
  if(model->opn + model->opn_un < 90.0)
    Lopn_cos_un2 = cos((model->opn + model->opn_un)/180.0*PI);
  else
    Lopn_cos_un2 = 0.0;

  /* generate clouds' radial location */
  if( flag_update == 1 )
  {
    pr = clouds_particles_perturb[which_particle_update];

    for(i=num_sh; i<parset.n_cloud_per_task; i++)
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial location of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      pr[i] = r;
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for(i=num_sh; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos_un2 + (Lopn_cos_un1-Lopn_cos_un2) * pow(gsl_rng_uniform(gsl_r), gam));
    //Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

    r = pr[i];
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then roate around z axis 
 * by an angle of Lphi
 */
  /*xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y - sin(Lthe)*cos(Lphi) * z;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb = sin(Lthe) * x + cos(Lthe) * z; */
    
    xb = cos_Lthe*cos_Lphi * x + sin_Lphi * y;
    yb =-cos_Lthe*sin_Lphi * x + cos_Lphi * y;
    zb = sin_Lthe * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;


    Vkep = sqrt(mbh/r);

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV = (gsl_ran_ugaussian(gsl_r) * sigthe_circ + 0.5) * PI;
      }
      else
      {
        if(fflow <= 0.5)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad + 1.0) *PI + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV = (gsl_ran_ugaussian(gsl_r) * sigthe_rad) * PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos_phi - Vph * sin_phi;
      vy = Vr * sin_phi + Vph * cos_phi;
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos_Lthe*cos_Lphi * vx + sin_Lphi * vy;
      vyb =-cos_Lthe*sin_Lphi * vx + cos_Lphi * vy;
      vzb = sin_Lthe * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
      vy = vyb;
      vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += gsl_ran_ugaussian(gsl_r) * sig_turb * Vkep;

      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = sqrt( (1.0 + V/C_Unit) / (1.0 - V/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;

      V += linecenter;

      tmp_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}

/*================================================================
 * model 8
 * binary model
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */

void transfun_1d_cloud_sample_model8(const void *pm, int flag_save)
{
  int i, flag_update = 0, num_sh;
  double *pr;
  BLRmodel8 *model = (BLRmodel8 *)pm;

  double L0, phi0, inc;
  L0 = exp(model->log_rLR) * (exp(model->log_Rin_1 + model->log_roi_1) + exp(model->log_Rin_2 + model->log_roi_2));
  phi0 = model->phi0 * PI / 180.0;
  inc = acos(model->cos_inc);

  /* BLR1*/
  double Rin, Rout, disp, cos_Theta_disk, gamma;
  Rin = exp(model->log_Rin_1);
  Rout = Rin * exp(model->log_roi_1);
  cos_Theta_disk = cos(model->Theta_disk_1 * PI / 180.0);
  gamma = model->gamma_1;
  disp = L0 * (1 - model->mu1);

  for (i = 0; i < num_params_radial_samp; i++)
  {
    if (which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the
  // initial step.
  if (force_update == 1 || which_parameter_update == -1)
  {
    flag_update = 1;
  }
  else
  {
    for (i = 0; i < num_params_radial_samp; i++)
    {
      if (which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  /* number of particles in first region */

  num_sh = (int)(parset.n_cloud_per_task * model->f1);

  /* generate clouds' radial location */
  if (flag_update == 1)
  {
    pr = clouds_particles_perturb[which_particle_update];
    for (i = 0; i < num_sh; i++)
    {
      pr[i] = Rin  + (Rout - Rin) * gsl_rng_uniform(gsl_r);
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  double r, phi, Lphi, Lthe, x, y, z, xb, yb, zb, dis;
  for (i = 0; i < num_sh; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
    Lphi = 2.0 * PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(cos_Theta_disk + (1.0 - cos_Theta_disk) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0 * PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi);
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe) * cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe) * sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    x = -cos(phi0) * (xb - disp) + sin(phi0) * yb;
    y = -sin(phi0) * (xb - disp) - cos(phi0) * yb;
    z = zb;

    dis = r + xb * sin(phi0) * sin(inc) + yb * cos(phi0) * sin(inc) - zb * cos(inc);
    tmp_tau[i] = dis;
    tmp_weight[i] = pow(r / Rin, gamma);

    if (flag_save && thistask == roottask)
    {
      if (i % (icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  /* BLR2 */
  Rin = exp(model->log_Rin_2);
  Rout = Rin * exp(model->log_roi_2);
  cos_Theta_disk = cos(model->Theta_disk_2 * PI / 180.0);
  gamma = model->gamma_2;
  disp = L0 * model->mu1;
  phi0 = PI + phi0;

  /* generate clouds' radial location */
  if (flag_update == 1)
  {
    pr = clouds_particles_perturb[which_particle_update];
    for (i = num_sh; i < parset.n_cloud_per_task; i++)
    {
      pr[i] = Rin + (Rout - Rin) * gsl_rng_uniform(gsl_r);
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for (i = 0; i < num_sh; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
    Lphi = 2.0 * PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(cos_Theta_disk + (1.0 - cos_Theta_disk) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0 * PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi);
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe) * cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe) * sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    x = -cos(phi0) * (xb - disp) + sin(phi0) * yb;
    y = -sin(phi0) * (xb - disp) - cos(phi0) * yb;
    z = zb;

    dis = r + xb * sin(phi0) * sin(inc) + yb * cos(phi0) * sin(inc) - zb * cos(inc);
    tmp_tau[i] = dis;
    tmp_weight[i] = pow(r / Rin, gamma);

    if (flag_save && thistask == roottask)
    {
      if (i % (icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }
  return;
}

/* 
 * This function caclulate 2d transfer function.
 */
void transfun_2d_cloud_sample_model8(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, flag_update = 0, num_sh, linecenter = 0.0;
  double *pr;
  double *pmodel = (double *)pm;
  BLRmodel8 *model = (BLRmodel8 *)pm;

  double L0, phi0, inc, Omega;
  L0 = exp(model->log_rLR) * (exp(model->log_Rin_1 + model->log_roi_1) + exp(model->log_Rin_2 + model->log_roi_2));
  phi0 = model->phi0 * PI / 180.0;
  inc = acos(model->cos_inc);
  Omega = sqrt(exp(model->log_mbh) / L0) / L0;

  /* BLR1*/
  double mbh, Rin, Rout, disp, cos_Theta_disk, gamma, Rs;
  mbh = exp(model->log_mbh) * model->mu1;
  Rin = exp(model->log_Rin_1);
  Rout = Rin * exp(model->log_roi_1);
  cos_Theta_disk = cos(model->Theta_disk_1 * PI / 180.0);
  gamma = model->gamma_1;
  disp = (1 - model->mu1) * L0;
  Rs = 3.0e11 * mbh / CM_PER_LD;

  if (parset.flag_linecenter != 0)
  {
    linecenter = pmodel[num_params_blr - num_params_linecenter - 1] * parset.linecenter_err;
  }

  if (force_update == 1 || which_parameter_update == -1)
  {
    flag_update = 1;
  }
  else
  {
    for (i = 0; i < num_params_radial_samp; i++)
    {
      if (which_parameter_update == params_radial_samp[i])
      {
        flag_update = 1;
        break;
      }
    }
  }

  num_sh = (int)(parset.n_cloud_per_task * model->f1);

  /* generate clouds' radial location */
  if (flag_update == 1)
  {
    pr = clouds_particles_perturb[which_particle_update];
    for (i = 0; i < num_sh; i++)
    {
      pr[i] = 1.5 * Rs + Rin + (Rout - Rin) * gsl_rng_uniform(gsl_r);
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  double r, phi, Lphi, Lthe, x, y, z, xb, yb, zb, dis;
  double Vkep, vx, vy, vz, vxb, vyb, vzb, V;
  double weight, g;
  for (i = 0; i < num_sh; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
    Lphi = 2.0 * PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(cos_Theta_disk + (1.0 - cos_Theta_disk) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0 * PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi);
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe) * cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe) * sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    x = -cos(phi0) * (xb - disp) + sin(phi0) * yb;
    y = -sin(phi0) * (xb - disp) - cos(phi0) * yb;
    z = zb;

    dis = r + xb * sin(phi0) * sin(inc) + yb * cos(phi0) * sin(inc) - zb * cos(inc);
    weight = pow(r / Rin, gamma);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    Vkep = sqrt(mbh / ( r - Rs));

    for (j = 0; j < parset.n_vel_per_cloud; j++)
    {
      vx = -Vkep * sin(phi);
      vy = Vkep * cos(phi);
      vz = 0.0;

      vxb = cos(Lthe) * cos(Lphi) * vx + sin(Lphi) * vy;
      vyb = -cos(Lthe) * sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      vx = -cos(phi0) * vxb + sin(phi0) * vyb - Omega * y;
      vy = -sin(phi0) * vxb - cos(phi0) * vyb + Omega * x;
      vz = vzb;

      V = -vy * sin(inc) - vz * cos(inc);
      if (fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999 * C_Unit * (V > 0.0 ? 1.0 : -1.0);

      g = sqrt((1.0 + V / C_Unit) / (1.0 - V / C_Unit)) / sqrt(1.0 - Rs / r); //relativistic effects
      V = (g - 1.0) * C_Unit;

      V += linecenter;
      tmp_vel[i * parset.n_vel_per_cloud + j] = V;

      if (flag_save && thistask == roottask)
      {
        if (i % (icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx * VelUnit, vy * VelUnit, vz * VelUnit, weight);
      }
    }
  }

  /* BLR2*/
  mbh = exp(model->log_mbh) * (1 - model->mu1);
  Rin = exp(model->log_Rin_2);
  Rout = Rin * exp(model->log_roi_2);
  cos_Theta_disk = cos(model->Theta_disk_2 * PI / 180.0);
  gamma = model->gamma_2;
  disp = model->mu1 * L0;
  Rs = 3.0e11 * mbh / CM_PER_LD;
  phi0 = PI + phi0;

  /* generate clouds' radial location */
  if (flag_update == 1)
  {
    pr = clouds_particles_perturb[which_particle_update];
    for (i = num_sh; i < parset.n_cloud_per_task; i++)
    {
      pr[i] = 1.5 * Rs + Rin + (Rout - Rin) * gsl_rng_uniform(gsl_r);
    }
  }
  else
  {
    pr = clouds_particles[which_particle_update];
  }

  for (i = num_sh; i < parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
    Lphi = 2.0 * PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(cos_Theta_disk + (1.0 - cos_Theta_disk) * gsl_rng_uniform(gsl_r));

    r = pr[i];
    phi = 2.0 * PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi);
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe) * cos(Lphi) * x + sin(Lphi) * y;
    yb = -cos(Lthe) * sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    x = -cos(phi0) * (xb - disp) + sin(phi0) * yb;
    y = -sin(phi0) * (xb - disp) - cos(phi0) * yb;
    z = zb;

    dis = r + xb * sin(phi0) * sin(inc) + yb * cos(phi0) * sin(inc) - zb * cos(inc);
    weight = pow(r / Rin, gamma);
    tmp_tau[i] = dis;
    tmp_weight[i] = weight;

    Vkep = sqrt(mbh / ( r - Rs));

    for (j = 0; j < parset.n_vel_per_cloud; j++)
    {
      vx = -Vkep * sin(phi);
      vy = Vkep * cos(phi);
      vz = 0.0;

      vxb = cos(Lthe) * cos(Lphi) * vx + sin(Lphi) * vy;
      vyb = -cos(Lthe) * sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      vx = -cos(phi0) * vxb + sin(phi0) * vyb - Omega * y;
      vy = -sin(phi0) * vxb - cos(phi0) * vyb + Omega * x;
      vz = vzb;

      V = -vy * sin(inc) - vz * cos(inc);
      if (fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999 * C_Unit * (V > 0.0 ? 1.0 : -1.0);

      g = sqrt((1.0 + V / C_Unit) / (1.0 - V / C_Unit)) / sqrt(1.0 - Rs / r); //relativistic effects
      V = (g - 1.0) * C_Unit;

      V += linecenter;
      tmp_vel[i * parset.n_vel_per_cloud + j] = V;

      if (flag_save && thistask == roottask)
      {
        if (i % (icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx * VelUnit, vy * VelUnit, vz * VelUnit, weight);
      }
    }
  }
  return;
}

void restart_action_1d(int iflag)
{
  FILE *fp;
  char str[200];
  int i, count;

  sprintf(str, "%s/data/clouds_%04d.txt", parset.file_dir, thistask);

  if(iflag == 0)  // write
    fp = fopen(str, "wb");
  else           // read
    fp = fopen(str, "rb");

  if(fp == NULL)
  {
    printf("# Cannot open file %s.\n", str);
  }

  if(iflag == 0)
  {
    printf("# Writing clouds at task %d.\n", thistask);
    for(i=0; i<parset.num_particles; i++)
    {
      count = fwrite(clouds_particles[i], sizeof(double), parset.n_cloud_per_task, fp);
      if(count < parset.n_cloud_per_task)
      {
        printf("# Error in writing clouds at task %d.\n", thistask);
      }
    }
  }
  else
  {
    printf("# Reading clouds at task %d.\n", thistask);
    for(i=0; i<parset.num_particles; i++)
    {
      count = fread(clouds_particles[i], sizeof(double), parset.n_cloud_per_task, fp);
      if(count < parset.n_cloud_per_task)
      {
        printf("# Error in reading clouds at task %d.\n", thistask);
      }
    }

  }
  fclose(fp);
}

void restart_action_2d(int iflag)
{
  restart_action_1d(iflag);
}

