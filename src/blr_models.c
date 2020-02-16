/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file blr_models.c
 *  \brief BLR models
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

/*================================================================
 * model 1
 * Brewer et al. (2011)'s model
 *
 * geometry: radial Gamma distribution
 * dynamics: elliptical orbits
 *================================================================
 */
/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void gen_cloud_sample_model1(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, L, E;
  double V, weight, rnd;
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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
      nc++;
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

    /* counter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif
    
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
      Vph = L/r; /* RM cannot distinguish the orientation of the rotation. */

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

      V = -vx;  //note the definition of the line-of-sight velocity. positive means a receding 
                // velocity relative to the observer.
      
      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
 * dynamics: elliptical orbits (Gaussian around circular orbits)
 *================================================================
 */
/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void gen_cloud_sample_model2(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe;
  double V, weight, rnd;
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

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum */    
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
      r = rin + sig * rnd;
      nc++;
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

    /* counter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif
    
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
      Vph = rhor * sin(rhotheta) * Vcirc; //RM cannot distinguish the orientation of the rotation.

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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

/*! 
 * This function calculate 2d transfer function at velocity grid "transv" and time grid "TransTau" . 
 * Note that time-lag grid is already set by parset.n_tau.
 */
void gen_cloud_sample_model3(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E;
  double V, weight, rnd;
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

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

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

// counter-rotate around y
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif
    
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
      Vph = L/r; //RM cannot distinguish the orientation of the rotation.

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
void gen_cloud_sample_model4(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E;
  double V, weight, rnd;
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

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum   */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));

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

    /* counter-rotate around y */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif
    
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
      Vph = sqrt(1.0-2.0*xi*xi) * L/r; //RM cannot distinguish the orientation of the rotation.

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
 * geometry: double power-law 
 * dynamics: elliptical orbits and inflow/outflow as in Pancoast's model
 *====================================================================
 */

/* 
 * This function caclulate 1d transfer function.
 */
void gen_cloud_sample_model5(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, V, Vr, Vph, Vkep, rhoV, theV;
  double cos_Lphi, sin_Lphi, cos_Lthe, sin_Lthe, cos_inc_cmp, sin_inc_cmp;
  double weight, rndr, rnd, rnd_frac, rnd_xi, frac1, frac2, ratio, Rs, g;
  double vx, vy, vz, vxb, vyb, vzb;
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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    cos_Lphi = cos(Lphi);
    sin_Lphi = sin(Lphi);
    cos_Lthe = cos(Lthe);
    sin_Lthe = sin(Lthe);
    
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
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
 * This function caclulate 2d transfer function.
 */
void gen_cloud_sample_model6(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double V, rhoV, theV, Vr, Vph, Vkep, Rs, g;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp;
  double weight, rnd, rnd_xi;
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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
      r = rin + sig * rnd;
      nc++;
    }
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

// counter-rotate around y, LOS is x-axis 
    /* x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
       y = yb;
       z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc); */
    
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
 * This function caclulate 2d transfer function.
 */
void gen_cloud_sample_model7(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc, num_sh;
  double r, phi, dis, Lopn_cos, cos_phi, sin_phi;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb, Rs, g, sig_turb;
  double V, rhoV, theV, Vr, Vph, Vkep;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp;
  double weight, rnd, rnd_xi;
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


  num_sh = (int)(parset.n_cloud_per_task * model->fsh);

  for(i=0; i<num_sh; i++)
  {
// generate a direction of the angular momentum of the orbit   
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * pow(gsl_rng_uniform(gsl_r), gam));
    sin_Lphi = sin(Lphi);
    cos_Lphi = cos(Lphi);
    sin_Lthe = sin(Lthe);
    cos_Lthe = cos(Lthe);

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
      r = rin + sig * rnd;
      nc++;
    }
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

// counter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
      r = rin + sig * rnd;
      nc++;
    }
    phi = 2.0*PI * gsl_rng_uniform(gsl_r);
    cos_phi = cos(phi);
    sin_phi = sin(phi);

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos_phi; 
    y = r * sin_phi;
    z = 0.0;

/* right-handed framework
 * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

// counter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    weight = 0.5 + k*(x/r);
    clouds_weight[i] = weight;
    
#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif

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

      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}


/* 
 * This function caclulate 2d transfer function.
 */
void gen_cloud_sample_model8(const void *pm, int flag_type, int flag_save)
{
  int i;
  double theta_min, theta_max, r_min, r_max, Rblr, Rv, mbh, alpha, gamma, xi, lambda, k;
  double rnd, r, r0, theta, phi, xb, yb, zb, x, y, z, l, weight, sin_inc_cmp, cos_inc_cmp, inc;
  double dis, density, vl, vesc, Vr, Vph, vx, vy, vz, vxb, vyb, vzb, V, rnd_xi, zb0, lmax, R;
  double v0=6.0/VelUnit;
  BLRmodel8 *model=(BLRmodel8 *)pm;

  theta_min = model->theta_min/180.0 * PI;
  theta_max = model->dtheta_max/180.0*PI + theta_min;
  theta_max = fmin(theta_max, 0.5*PI);

  r_min = exp(model->r_min);
  r_max = model->fr_max * r_min;
  if(r_max > 0.5*rcloud_max_set) r_max = 0.5*rcloud_max_set;
  
  gamma = model->gamma;
  alpha = model->alpha;
  lambda = model->lamda;
  
  k = model->k;
  xi = model->xi;

  Rv = exp(model->Rv);
  Rblr = exp(model->Rblr);
  if(Rblr < r_max) Rblr = r_max;
  inc = acos(model->inc);
  
  mbh = exp(model->mbh);

  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    rnd = gsl_rng_uniform(gsl_r);
    r0 = r_min +  rnd * (r_max - r_min);
    theta = theta_min + (theta_max - theta_min) * pow(rnd, gamma); 
    /* the maximum allowed value of l */
    lmax = -r0*sin(theta) + sqrt(r0*r0*sin(theta)*sin(theta) + (Rblr*Rblr - r0*r0));
    l = gsl_rng_uniform(gsl_r) * lmax;

    r = l * sin(theta) + r0;
    if(gsl_rng_uniform(gsl_r) < 0.5)
      zb = l * cos(theta);
    else
      zb =-l * cos(theta);
      
    phi = gsl_rng_uniform(gsl_r) * 2.0*PI;

    xb = r * cos(phi);
    yb = r * sin(phi);

    zb0 = zb;
    rnd_xi = gsl_rng_uniform(gsl_r);
    if( (rnd_xi < 1.0 - xi) && zb0 < 0.0)
      zb = -zb;

    vesc = sqrt(2.0*mbh/r0);
    vl = v0 + (vesc - v0) * pow(l/Rv, alpha)/(1.0 + pow(l/Rv, alpha));
    density = pow(r0, lambda)/vl;

    // counter-rotate around y, LOS is x-axis 
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    R = sqrt(r*r + zb*zb);
    weight = 0.5 + k*(x/R);
    clouds_weight[i] = weight * density;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = R - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = R - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = R - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = R - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif

    Vr = vl * sin(theta);
    vzb = vl * cos(theta);
    Vph = sqrt(mbh*r0) / r;

    vxb = Vr * cos(phi) - Vph * sin(phi);
    vyb = Vr * sin(phi) + Vph * cos(phi);

    if(zb < 0.0)
      vzb = -vzb;

    vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
    vy = vyb;
    vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

    V = -vx;  //note the definition of the line-of-sight velocity. postive means a receding 
                // velocity relative to the observer.
    clouds_vel[i] = V;

    if(flag_save && thistask==roottask)
    {
      if(i%(icr_cloud_save) == 0)
        fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
    }
  }

  return;
}

void gen_cloud_sample_model9(const void *pm, int flag_type, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn, Lopn_cos, u;
  double x, y, z, xb, yb, zb, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, a, s, rin, sig;
  double Lphi, Lthe, Vkep, Rs, g;
  double V, weight, rnd, sin_inc_cmp, cos_inc_cmp;
  BLRmodel9 *model = (BLRmodel9 *)pm;
  double Vr, Vph, mbh;
  
  //Lopn_cos = cos(model->opn*PI/180.0);
  Lopn = model->opn*PI/180.0;
  inc = acos(model->inc);
  beta = model->beta;
  F = model->F;
  mu = exp(model->mu);

  mbh = exp(model->mbh);
  Rs = 3.0e11*mbh/CM_PER_LD; // Schwarzchild radius in a unit of light-days

  a = 1.0/beta/beta;
  s = mu/a;
  rin = mu*F + Rs;
  sig = (1.0-F)*s;

  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum     
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    //Lthe = acos(Lopn_cos + (1.0-Lopn_cos) * gsl_rng_uniform(gsl_r));
    Lthe = gsl_rng_uniform(gsl_r) * Lopn;

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
//    r = mu * F + (1.0-F) * gsl_ran_gamma(gsl_r, 1.0/beta/beta, beta*beta*mu);
        r = rin + sig * rnd;
      nc++;
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

    /* counter-rotate around y */
    //x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    //y = yb;
    //z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);
    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    weight = 1.0;
    clouds_weight[i] = weight;

#ifndef SA
    dis = r - x;
    clouds_tau[i] = dis;
    if(flag_type == 1)
    {
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
    }
#else
  switch(flag_type) 
  {
    case 1:   /* 1D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\n", x, y, z);
      }
      continue;
      break;
    
    case 2:  /* 2D RM */
      dis = r - x;
      clouds_tau[i] = dis;
      break;

    case 3: /* SA */
      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 4: /* 1D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
    
    case 5: /* 2D RM + SA */
      dis = r - x;
      clouds_tau[i] = dis;

      clouds_alpha[i] = y;
      clouds_beta[i] = z;
      break;
  }
  
#endif
    
    /* velocity  
     * note that a cloud moves in its orbit plane, whose direction
     * is determined by the direction of its angular momentum.
     */
    Vkep = sqrt(mbh/r);
      
    for(j=0; j<parset.n_vel_per_cloud; j++)
    {
      Vr = 0.0;
      Vph = Vkep; /* RM cannot distinguish the orientation of the rotation. */

      vx = Vr * cos(phi) - Vph * sin(phi);
      vy = Vr * sin(phi) + Vph * cos(phi);
      vz = 0.0;     

    /*vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy - sin(Lthe)*cos(Lphi) * vz;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy + sin(Lthe)*sin(Lphi) * vz;
      vzb = sin(Lthe) * vx + cos(Lthe) * vz;*/

      vxb = cos(Lthe)*cos(Lphi) * vx + sin(Lphi) * vy;
      vyb =-cos(Lthe)*sin(Lphi) * vx + cos(Lphi) * vy;
      vzb = sin(Lthe) * vx;
    
      //vx = vxb * cos(PI/2.0-inc) + vzb * sin(PI/2.0-inc);
      //vy = vyb;
      //vz =-vxb * sin(PI/2.0-inc) + vzb * cos(PI/2.0-inc);
      vx = vxb * cos_inc_cmp + vzb * sin_inc_cmp;
      vy = vyb;
      vz =-vxb * sin_inc_cmp + vzb * cos_inc_cmp;

      V = -vx;  //note the definition of the line-of-sight velocity. positive means a receding 
                // velocity relative to the observer.
      
      if(fabs(V) >= C_Unit) // make sure that the velocity is smaller than speed of light
        V = 0.9999*C_Unit * (V>0.0?1.0:-1.0);

      g = (1.0 + V/C_Unit) / sqrt( (1.0 - V*V/C_Unit/C_Unit) ) / sqrt(1.0 - Rs/r); //relativistic effects
      V = (g-1.0)*C_Unit;
      
      clouds_vel[i*parset.n_vel_per_cloud + j] = V;

      if(flag_save && thistask==roottask)
      {
        if(i%(icr_cloud_save) == 0)
          fprintf(fcloud_out, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", x, y, z, vx*VelUnit, vy*VelUnit, vz*VelUnit, weight);
      }
    }
  }

  return;
}

void restart_action_1d(int iflag)
{
  /* at present, nothing needs to do */
  /*
  FILE *fp;
  char str[200];

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
    printf("# Writing restart at task %d.\n", thistask);
  }
  else
  {
    printf("# Reading restart at task %d.\n", thistask);
  }
  fclose(fp);*/
  return;
}

void restart_action_2d(int iflag)
{
  restart_action_1d(iflag);
}

#ifdef SA
void restart_action_sa(int iflag)
{
  restart_action_1d(iflag);
}
#endif