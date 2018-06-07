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
 * Note that the light curves has been obtained in advance 
 */
void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl)
{
  int i, j, k;
  double fline, fcon, tl, tc, tau, A, Ag, ftrend;
  double *pmodel = (double *)pm;

  A=exp(pmodel[0]);
  Ag=pmodel[1];

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
        //fcon = mean; 
        fcon = con_q[0];
        for(k=1; k < nq; k++)/*  beyond the range, set to be the long-term trend */
        {
          fcon += con_q[k] * pow(tc, k);
        }
      }
      
      // add different trend in continuum and emission
      if(parset.flag_trend_diff)
      {
        ftrend = pmodel[num_params_blr + 4 + parset.flag_trend] * (tc - 0.5*(Tcon_data[0] + Tcon_data[n_con_data-1]));
        fcon += ftrend;
      }

      if(fcon > 0.0)
      {
        fline += Trans1D[j] * fcon * pow(fcon, Ag);     /*  line response */
      }	
    }
    fline *= dTransTau * A;
    Fl[i] = fline;
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
  double tau, tl, tc, fcon, A, Ag, ftrend, fnarrow;
  double *pmodel = (double *)pm;

  A=exp(pmodel[0]);
  Ag=pmodel[1];

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
        //fcon = mean; 
        fcon = con_q[0];
        for(m=1; m < nq; m++)/*  beyond the range, set to be the long-term trend */
        {
            fcon += con_q[m] * pow(tc, m);
        }
      }

      // add different trend in continuum and emission line
      if(parset.flag_trend_diff)
      {
        ftrend = pmodel[num_params_blr + 4 + parset.flag_trend] * (tc - 0.5*(Tcon_data[0] + Tcon_data[n_con_data-1]));
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

  // add intrinsic narrow line
  if(parset.flag_narrowline != 0)
  {
    double flux, width, shift;
    if(parset.flag_narrowline == 1)
    {
      flux = parset.flux_narrowline;
      width = parset.width_narrowline;
      shift = parset.shift_narrowline;
    }
    else if(parset.flag_narrowline == 2)
    {
      flux = parset.flux_narrowline + pmodel[num_params_blr-num_params_res-1-3] * parset.flux_narrowline_err;
      width = parset.width_narrowline + pmodel[num_params_blr-num_params_res-1-2] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr-num_params_res-1-1] * parset.shift_narrowline_err;
    }
    else
    {
      flux = exp(pmodel[num_params_blr-num_params_res-1-3]);
      width = parset.width_narrowline + pmodel[num_params_blr-num_params_res-1-2] * parset.width_narrowline_err;
      shift = parset.shift_narrowline + pmodel[num_params_blr-num_params_res-1-1] * parset.shift_narrowline_err;
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

// smooth the line profile
  line_gaussian_smooth_2D_FFT(transv, fl2d, nl, nv, pm);
}

/*================================================================
 * model 1
 * Brewer et al. (2011)'s model
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_direct_model1(const void *pm, int flag_save)
{
  FILE *fcloud_out;
  int i, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe;
  double Anorm, weight, rnd;
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
  
  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

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

  /* reset transfer function */
  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //	continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }
  
  if(flag_save && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}


/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_direct_model1(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idV, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, L, E, vcloud_max, vcloud_min;
  double dV, V, Anorm, weight, rnd;
  BLRmodel1 *model = (BLRmodel1 *)pm;
  FILE *fcloud_out;
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
  

  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    if( flag_update == 1 )
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

    //if(dis< parset.tau_min_set || dis>= parset.tau_max_set + dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    
// velocity  
// note that a cloud moves in its orbit plane, whose direction
// is determined by the direction of its angular momentum.
    Emin = - mbh/r;
    //Ecirc = 0.5 * Emin;
    //Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      chi = lambda * gsl_ran_ugaussian(gsl_r);
      E = Emin / (1.0 + exp(-chi));

      Lmax = sqrt(2.0 * r*r * (E + mbh/r));      

      if(lambda>1.0e-2)   //make sure that exp is caculatable.
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

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}


/*================================================================
 * model 2
 *================================================================
 */
/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_direct_model2(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idV, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, vcloud_max, vcloud_min;
  double dV, V, Anorm, weight, rnd;
  BLRmodel2 *model = (BLRmodel2 *)pm;
  FILE *fcloud_out;
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

  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    // which_parameter_update == 1 ==> beta is being updated.
    if( flag_update == 1 )
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

    //if(dis< parset.tau_min_set || dis>= parset.tau_max_set + dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    
// velocity  
// note that a cloud moves in its orbit plane, whose direction
// is determined by the direction of its angular momentum.
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

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}


/*====================================================================
 * model 3
 *====================================================================
 */

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_direct_model3(const void *pm, int flag_save)
{
  FILE *fcloud_out;
  int i, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe;
  double Anorm, weight, rnd;
  BLRmodel3 *model = (BLRmodel3 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  alpha = model->alpha;         
  F = exp(model->F);
  Rin = exp(model->Rin);                 /* mean radius */
  k = model->k;  
  
  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  /* reset transfer function */
  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));
    
    // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the 
    // initial step.
    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }
  
  if(flag_save && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}

/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_direct_model3(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idV, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, vcloud_max, vcloud_min;
  double dV, V, Anorm, weight, rnd;
  BLRmodel3 *model = (BLRmodel3 *)pm;
  FILE *fcloud_out;
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
  
  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    // which_parameter_update == 1 ==> beta is being updated.
    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

    //if(dis< parset.tau_min_set || dis>= parset.tau_max_set + dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    
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
      u = gsl_rng_uniform(gsl_r);
      Vr = Vr * (u<q?-1.0:1.0);
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

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}

/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void transfun_2d_cloud_direct_model4(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idV, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, vcloud_max, vcloud_min;
  double dV, V, Anorm, weight, rnd;
  BLRmodel4 *model = (BLRmodel4 *)pm;
  FILE *fcloud_out;
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
  
  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }
  
  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

    // which_parameter_update == 1 ==> beta is being updated.
    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        if(fabs(1.0-alpha) > 1.0e-4)
           rnd = pow( gsl_rng_uniform(gsl_r) * ( pow(F, 1.0-alpha) - 1.0) + 1.0,  1.0/(1.0-alpha) );
        else
           rnd = exp( gsl_rng_uniform(gsl_r) * log(F) );
        r = Rin * rnd ;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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

// conter-rotate around y
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis< parset.tau_min_set || dis>= parset.tau_max_set + dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    
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
      u = gsl_rng_uniform(gsl_r);
      Vr = Vr * (u<q?-1.0:1.0);
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

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;
 
      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}

/*====================================================================
 * model 5
 *====================================================================
 */

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_direct_model5(const void *pm, int flag_save)
{
  FILE *fcloud_out;
  int i, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double Lphi, Lthe;
  double Anorm, weight, rndr, rnd, rnd_xi, rnd_frac, frac1, frac2, ratio;
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
  

  frac1 = 1.0/(alpha+1.0) * (1.0 - pow(Fin, alpha+1.0));
  frac2 = 1.0/(alpha-1.0) * (1.0 - pow(Fout, -alpha+1.0));
  ratio = frac1/(frac1 + frac2);

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }

  /* reset transfer function */
  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    
    // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the 
    // initial step.
    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
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
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;

    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }
  
  if(flag_save && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}

/* 
 * This function caclulate 1d transfer function.
 */
void transfun_2d_cloud_direct_model5(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  FILE *fcloud_out;
  int i, j, idt, idV, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, rf;
  double Lphi, Lthe, V, Vr, Vph, Vkep, rhoV, theV;
  double Anorm, weight, rndr, rnd, rnd_frac, rnd_xi,rnd_flow, frac1, frac2, ratio, fe;
  double vx, vy, vz, vxb, vyb, vzb, dV, vcloud_max, vcloud_min;
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
  rf = exp(model->rf);
  sigr_circ = exp(model->sigr_circ);
  sigthe_circ = exp(model->sigthe_circ);
  sigr_rad = exp(model->sigr_rad);
  sigthe_rad = exp(model->sigthe_rad);
  theta_rot = model->theta_rot*PI/180.0;

  frac1 = 1.0/(alpha+1.0) * (1.0 - pow(Fin,   alpha+1.0));
  frac2 = 1.0/(alpha-1.0) * (1.0 - pow(Fout, -alpha+1.0));
  ratio = frac1/(frac1 + frac2);

  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

  for(i=0; i<num_params_radial_samp;i++)
  {
    if(which_parameter_update == params_radial_samp[i])
    {
      flag_update = 1;
      break;
    }
  }
  if(force_update == 1 || which_parameter_update == -1 )
  {
    flag_update = 1;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    
    // "which_parameter_update = -1" means that all parameters are updated, usually occurs at the 
    // initial step.
    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
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
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// counter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;

    Vkep = sqrt(mbh/r);
     
    if(r >= rf)
      fe = fellip;
    else
      fe = 0.0;

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);
      rnd_flow = gsl_rng_uniform(gsl_r);

      if(rnd < fe)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV =  gsl_ran_ugaussian(gsl_r) * sigthe_circ + PI/2.0;
      }
      else
      {
        if(rnd_flow < fflow)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}


/*================================================================
 * model 6
 * Pancoast et al. (2014)'s model
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_direct_model6(const void *pm, int flag_save)
{
  FILE *fcloud_out;
  int i, idt, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
  double Lphi, Lthe;
  double Anorm, weight, rnd, rnd_xi;
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
  
  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

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

  /* reset transfer function */
  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    
    xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }
  
  if(flag_save && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}

/* 
 * This function caclulate 2d transfer function.
 */
void transfun_2d_cloud_direct_model6(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  FILE *fcloud_out;
  int i, j, idt, idV, nc, flag_update=0;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb, vcloud_min, vcloud_max;
  double V, dV, rhoV, theV, Vr, Vph, Vkep, Rs, g, Vt;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe;
  double Anorm, weight, rnd, rnd_xi, rnd_flow;
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
  rin=mu*F;
  sig=(1.0-F)*s;
  
  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }

    r += Rs; //add Schwarzschild radius
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
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);

    Vkep = sqrt(mbh/r);
    
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);
      rnd_flow = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV =  gsl_ran_ugaussian(gsl_r) * sigthe_circ + PI/2.0;
      }
      else
      {
        if(rnd_flow < fflow)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;    
      
    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;


      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.

      V += gsl_ran_ugaussian(gsl_r) * sig_turb * Vkep; // add turbulence velocity

      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = sqrt( (1.0 + V/C_Unit) / (1.0 - V/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;
      
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}

/*================================================================
 * model 7
 * shadowed model
 *================================================================
 */
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_direct_model7(const void *pm, int flag_save)
{
  FILE *fcloud_out;
  int i, idt, nc, flag_update=0, num_sh;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
  double Lphi, Lthe;
  double Anorm, weight, rnd, rnd_xi;
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
  
  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

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

  /* reset transfer function */
  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  // number of particles in first region
  
  num_sh = (int)(parset.n_cloud_per_task * model->fsh);

  for(i=0; i<num_sh; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    
    xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
  }

  //second BLR region
  rin = exp(model->mu_un) * model->F_un;
  a = 1.0/(model->beta_un * model->beta_un);
  sig = exp(model->mu_un) * (1.0-model->F_un)/a;
  double Lopn_cos_un1, Lopn_cos_un2;
  Lopn_cos_un1 = Lopn_cos;
  if(model->opn + model->opn_un < 90.0)
    Lopn_cos_un2 = cos((model->opn + model->opn_un)/180.0*PI);
  else
    Lopn_cos_un2 = 0.0;

  for(i=num_sh; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos_un2 + (Lopn_cos_un1-Lopn_cos_un2) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    
    xb = cos(Lthe)*cos(Lphi) * x + sin(Lphi) * y;
    yb =-cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y;
    zb = sin(Lthe) * x;

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;
    Trans1D[idt] += weight;

    if(flag_save && thistask==roottask)
    {
      fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
    }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      Trans1D[i] = 0.0;
    }
  }
  
  if(flag_save && thistask==roottask)
  {
    fclose(fcloud_out);
  }

  return;
}

/* 
 * This function caclulate 2d transfer function.
 */
void transfun_2d_cloud_direct_model7(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  FILE *fcloud_out;
  int i, j, idt, idV, nc, flag_update=0, num_sh;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb, vcloud_min, vcloud_max, fe;
  double V, dV, rhoV, theV, Vr, Vph, Vkep;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot;
  double Lphi, Lthe;
  double Anorm, weight, rnd, rnd_xi, rnd_flow;
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

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;
  
  dV =(transv[1] - transv[0]); // velocity grid width

  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]=0.0;   // cleanup of transfer function

  vcloud_max = -DBL_MAX;
  vcloud_min = DBL_MAX;

  if(flag_save && thistask == roottask)
  {
    char fname[200];
    sprintf(fname, "%s/%s", parset.file_dir, parset.cloud_out_file);
    fcloud_out = fopen(fname, "w");
    if(fcloud_out == NULL)
    {
      fprintf(stderr, "# Error: Cannot open file %s\n", fname);
      exit(-1);
    }
  }

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
  
  num_sh = (int)(parset.n_cloud_per_task * model->fsh);
  for(i=0; i<num_sh; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;


    Vkep = sqrt(mbh/r);

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);
      rnd_flow = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV =  gsl_ran_ugaussian(gsl_r) * sigthe_circ + PI/2.0;
      }
      else
      {
        if(rnd_flow < fflow)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
    }
  }

  //second BLR region
  fellip = model->fellip_un;
  fflow = model->fflow_un;
  rin = exp(model->mu_un) * model->F_un;
  a = 1.0/(model->beta_un * model->beta_un);
  sig = exp(model->mu_un) * (1.0-model->F_un)/a;
  double Lopn_cos_un1, Lopn_cos_un2;
  Lopn_cos_un1 = Lopn_cos;
  if(model->opn + model->opn_un < 90.0)
    Lopn_cos_un2 = cos((model->opn + model->opn_un)/180.0*PI);
  else
    Lopn_cos_un2 = 0.0;

  for(i=num_sh; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos_un2 + (Lopn_cos_un1-Lopn_cos_un2) * pow(gsl_rng_uniform(gsl_r), gam));
    //Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));

    if( flag_update == 1 ) //
    {
      nc = 0;
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        if(nc > 1000)
        {
          printf("# Error, too many tries in generating ridial localtion of clouds.\n");
          exit(0);
        }
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//      r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
        nc++;
      }
      clouds_particles_perturb[which_particle_update][i] = r;
    }
    else
    {
      r = clouds_particles[which_particle_update][i];
    }
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
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

// conter-rotate around y, LOS is x-axis 
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;

    //if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    //  continue;
    idt = (dis - parset.tau_min_set)/dTransTau;
    if(idt < 0)
      idt = 0;
    if(idt >= parset.n_tau)
      idt = parset.n_tau - 1;

    weight = 0.5 + k*(x/r);
    //weight = 0.5 + k * x/sqrt(x*x+y*y);
    //Trans1D[idt] += pow(1.0/r, 2.0*(1 + gam)) * weight;


    Vkep = sqrt(mbh/r);

    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      rnd = gsl_rng_uniform(gsl_r);
      rnd_flow = gsl_rng_uniform(gsl_r);

      if(rnd < fellip)
      {
        rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_circ  + 1.0) * Vkep;
        theV =  gsl_ran_ugaussian(gsl_r) * sigthe_circ + PI/2.0;
      }
      else
      {
        if(rnd_flow < fflow)
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + theta_rot;
        }
        else
        {
          rhoV = (gsl_ran_ugaussian(gsl_r) * sigr_rad  + 1.0) * Vkep;
          theV =  gsl_ran_ugaussian(gsl_r) * sigthe_rad + PI + theta_rot;
        }
      }
      
      Vr = sqrt(2.0) * rhoV * cos(theV);
      Vph = rhoV * sin(theV);

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;

      if((rnd_xi < 1.0-xi) && zb0 < 0.0)
        vzb = -vzb;
    
      vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      vy = vyb;
      vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);

      vcloud_max = fmax(vx, vcloud_max);
      vcloud_min = fmin(vx, vcloud_min);

      V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;

      idV = (V - transv[0])/dV;

      //trans2d[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + gam)) * weight;
      trans2d[idt*n_vel + idV] += weight;

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }
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
    printf(" Warning, zero transfer function at task %d.\n", thistask);
    for(i=0; i<parset.n_tau; i++)
    {
      for(j=0; j<n_vel; j++)
      {
        trans2d[i*n_vel+j] = 0.0;
      }
    }
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}


void restart_clouds_1d(int iflag)
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

void restart_clouds_2d(int iflag)
{
  restart_clouds_1d(iflag);
}

