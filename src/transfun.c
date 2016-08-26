/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
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

/* note that the light curves has been obtained in advance */

void calculate_line_from_blrmodel(const void *pm, double *Tl, double *Fl, int nl)
{
  int i, j;
  double fline, fcon, tl, tc, tau, A;
  BLRmodel *model = (BLRmodel *)pm;

  A=exp(model->A);
  
  //dtau = (TransTau[parset.n_tau-1] - TransTau[0])/(parset.n_tau-1);

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
  		  fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, tc, gsl_acc);
  		//	fline += Trans1D[j] * fcon * pow(fabs(fcon), pmodel.Ag);
        fline += Trans1D[j] * fcon;
  	  }
  	}
  	fline *= dTransTau * A;
  	Fl[i] = fline;
  }

}

void transfun_1d_cloud_direct(const void *pm)
{
  int i, idx;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb;
  double inc, F, beta, mu, k, a, s;
  double Lphi, Lthe;
  double Anorm, weight, rnd;
  BLRmodel *model = (BLRmodel *)pm;

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = model->inc * PI/180.0;
  beta = model->beta;
  F = model->F;
  mu = exp(model->mu);
  k = model->k;  

  a = 1.0/beta/beta;
  s = mu/a;
  

  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] = 0.0;
  }
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    if(Lopn_cos<1.0)
    {
      Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));
    }
    else
      Lthe = 0.0;
    
    // when the lastest MH move is not accepted, update clouds
    if(beta !=beta_old_particles[which_particle_update] || which_parameter_update == -1) //
    {
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//        r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, mu/beta, beta);
//        r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = mu * F + (1.0-F) * s * rnd;
      }
      clouds_particles[which_particle_update][i] = rnd;
    }
    else
    {
      rnd = clouds_particles[which_particle_update][i];
      r = mu * F + (1.0-F) * s * rnd;
    }
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe)*cos(Lphi) * x - sin(Lphi) * y + sin(Lthe)*cos(Lphi) * z;
    yb = cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =-sin(Lthe) * x + cos(Lthe) * z;

// rotate around y
    x = xb * cos(inc) - zb * sin(inc);
    y = yb;
    z = xb * sin(inc) + zb * cos(inc);

    dis = r - z;

    if(dis<parset.tau_min_set || dis>=parset.tau_max_set+dTransTau)
    	continue;
    //weight = 0.5 + k*(z/r);
    weight = 0.5 + k * x/sqrt(x*x+y*y);
    idx = (dis - parset.tau_min_set)/dTransTau;
    //Trans1D[idx] += pow(1.0/r, 2.0*(1 + model.Ag)) * weight;
    Trans1D[idx] += weight;
  }

  // record the previous beta to save comupational time
  if(which_parameter_update == 1 || which_parameter_update == -1)
    beta_old_particles[which_particle_update] = beta;

  Anorm = 0.0;
  for(i=0;i<parset.n_tau;i++)
  {
    Anorm += Trans1D[i] * dTransTau;
  }

  for(i=0; i<parset.n_tau; i++)
  {
    Trans1D[i] /= Anorm;
  }

  return;
}

void calculate_line2d_from_blrmodel(const void *pm, const double *Tl, const double *transv, const double *trans2d, 
                                              double *fl2d, int nl, int nv)
{
  int i, j, k;
  double fline, tau, tl, tc, fcon, A;
  BLRmodel *model = (BLRmodel *)pm;

  A=exp(model->A);

  for(j=0;j<nl; j++)
  {
    tl = Tl[j];
    for(i=0; i<nv; i++)
    {
      fline = 0.0;
      for(k=0; k<parset.n_tau; k++)
      {
        tau = TransTau[k];
        tc = tl - tau;
        if(tc>=Tcon_min && tc <=Tcon_max)
        {
          fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, tc, gsl_acc);
        //  fline += Trans2D[k*n_vel+i] * fcon * pow(fabs(fcon), pmodel.Ag);
          fline += trans2d[k*nv+i] * fcon;
        }
      }
      fline *= dTransTau * A ;
      fl2d[j*nv + i] = fline;
    }
  }

// smooth out the line profile
  line_gaussian_smooth_2D_FFT(transv, fl2d, nl, nv);
}

/* time-lag grid is already set by parset.n_tau 
 * calculate transfer function at velocity grid "transv" and time grid "TransTau" 
 */
void transfun_2d_cloud_direct(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, idV, idt;
  double vrange, r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s;
  double Lphi, Lthe, L, E, vcloud_max, vcloud_min;
  double dV, V, Anorm, weight, rnd;
  BLRmodel *model = (BLRmodel *)pm;
  FILE *fcloud_out;
  double Emin, Lmax, Vr, Vr2, Vph, mbh, chi, lambda, q;
  

  Lopn_cos = cos(model->opn*PI/180.0);
  inc = model->inc * PI/180.0;
  beta = model->beta;
  F = model->F;
  mu = exp(model->mu);
  k = model->k;

  a = 1.0/beta/beta;
  s = mu/a;

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

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    if(Lopn_cos<1.0)
    {
      Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));
    }
    else
      Lthe = 0.0;

    if(beta !=beta_old_particles[which_particle_update]) // beta updated
    {
      r = rcloud_max_set+1.0;
      while(r>rcloud_max_set || r<rcloud_min_set)
      {
        rnd = gsl_ran_gamma(gsl_r, a, 1.0);
//        r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, mu/beta, beta);
//        r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = mu * F + (1.0-F) * s * rnd;
      } 
      clouds_particles[which_particle_update][i] = rnd;
    }
    else
    {
      rnd = clouds_particles[which_particle_update][i];
      r = mu * F + (1.0-F) * s * rnd;
    }
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);

    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

    xb = cos(Lthe)*cos(Lphi) * x - sin(Lphi) * y + sin(Lthe)*cos(Lphi) * z;
    yb = cos(Lthe)*sin(Lphi) * x + cos(Lphi) * y + sin(Lthe)*sin(Lphi) * z;
    zb =-sin(Lthe) * x + cos(Lthe) * z;

// rotate around y
    x = xb * cos(inc) - zb * sin(inc);
    y = yb;
    z = xb * sin(inc) + zb * cos(inc);

    dis = r - z;

    if(dis< parset.tau_min_set || dis>= parset.tau_max_set + dTransTau)
      continue;
    //weight = 0.5 + k*(z/r);
    weight = 0.5 + k * x/sqrt(x*x+y*y);
    idt = (dis - parset.tau_min_set)/dTransTau;

// velocity  
// note that a cloud moves in its orbit plane, whose direction
// is determined by the direction of its angular momentum.
    Emin = - mbh/r;
    //Ecirc = 0.5 * Emin;
    //Lcirc = sqrt(2.0 * r*r * (Ecirc + mbh/r));
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      chi = lambda * gsl_ran_gaussian(gsl_r, 1.0);
      E = Emin / (1.0 + exp(-chi));
      Lmax = sqrt(2.0 * r*r * (E + mbh/r));
      if(lambda>1.0e-2)   //make sure that exp is caculatable.
        L = Lmax * lambda * log( (exp(1.0/lambda) - 1.0) * gsl_rng_uniform(gsl_r) + 1.0 );
      else
        L = Lmax + lambda * log(gsl_rng_uniform(gsl_r));
 
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
      u = gsl_rng_uniform(gsl_r);
      Vph = L/r * (u<0.5?-1.0:1.0);

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

      vxb = cos(Lthe)*cos(Lphi) * vx - sin(Lphi) * vy + sin(Lthe)*cos(Lphi) * vz;
      vyb = cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb =-sin(Lthe) * vx + cos(Lthe) * vz;
    
      vx = vxb * cos(inc) - vzb * sin(inc);
      vy = vyb;
      vz = vxb * sin(inc) + vzb * cos(inc);

      vcloud_max = fmax(vz, vcloud_max);
      vcloud_min = fmin(vz, vcloud_min);

      if(flag_save && thistask==roottask)
      {
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx, vy, vz);
      }

      V = -vz;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
      if(V<transv[0] || V>=transv[n_vel-1]+dV)
        continue;
      idV = (V - transv[0])/dV;
      //Trans2D[idt*n_vel + idV] += pow(1.0/r, 2.0*(1 + model.Ag)) * weight;
      trans2d[idt*n_vel + idV] += weight;
    }
  }
  
  // record the previous beta
  if(which_parameter_update == 1 || which_parameter_update == -1)
    beta_old_particles[which_particle_update] = beta;

  Anorm = 0.0;
  for(i=0; i<parset.n_tau; i++)
    for(j=0; j<n_vel; j++)
    {
      Anorm += trans2d[i*n_vel+j] * dV * dTransTau;
    }

  for(i=0; i<parset.n_tau; i++)
  {
    for(j=0; j<n_vel; j++)
    {
      trans2d[i*n_vel+j] /= Anorm;
    }
  }

  vrange = vcloud_max - vcloud_min;
  if(vrange<=0.0)
  {
    fprintf(stderr, "No velocity-resolved lag for given BLR model!\n");
    //printf("%f %f %f %f %f\n", model->inc, model->opn, model->lambda, model->F, model->beta);
    //printf("%f %f %f\n", model->k, model->q, model->mu);
    //printf("%e %e %f %f %f\n", vcloud_max, vcloud_min, V, Vr, Vph);
  }

  if(flag_save && thistask == roottask)
    fclose(fcloud_out);

  return;
}