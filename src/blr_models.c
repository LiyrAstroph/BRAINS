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
/* 
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model1(const void *pm, int flag_save)
{
  int i, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe;
  double weight, rnd;
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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */
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

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

   /* right-handed framework
    * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

    /* counter-rotate around y, LOS is x-axis */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
  int i, j, nc;
  double r, phi, dis, Lopn_cos, u;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, L, E, linecenter = 0.0;
  double V, weight, rnd;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;
    
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
      
      V += linecenter;
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
void transfun_2d_cloud_sample_model2(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, beta, mu, k, a, s, rin, sig;
  double Lphi, Lthe, linecenter = 0.0;
  double V, weight, rnd;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;
    
    
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

      V += linecenter;
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

/* 
 * This function calculate 1d transfer function.
 */
void transfun_1d_cloud_sample_model3(const void *pm, int flag_save)
{
  int i, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe;
  double weight, rnd;
  BLRmodel3 *model = (BLRmodel3 *)pm;

  Lopn_cos = cos(model->opn*PI/180.0); /* cosine of openning angle */
  inc = acos(model->inc);         /* inclination angle in rad */
  alpha = model->alpha;         
  F = exp(model->F);
  Rin = exp(model->Rin);                 /* mean radius */
  k = model->k;  
  
  if(F*Rin > rcloud_max_set)
    F = rcloud_max_set/Rin;

  for(i=0; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit */
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

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
    z = 0.0;

    /* right-handed framework
     * first rotate around y axis by an angle of Lthe, then rotate around z axis 
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

    /* counter-rotate around y, LOS is x-axis */
    x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc);

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, linecenter=0.0;
  double V, weight, rnd;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;
    
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

      V += linecenter;
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
void transfun_2d_cloud_sample_model4(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double inc, F, alpha, Rin, k;
  double Lphi, Lthe, L, E, linecenter=0.0;
  double V, weight, rnd;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;
    
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

      V += linecenter;
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
void transfun_1d_cloud_sample_model5(const void *pm, int flag_save)
{
  int i, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double Lphi, Lthe, cos_Lphi, sin_Lphi, cos_Lthe, sin_Lthe, cos_inc_cmp, sin_inc_cmp;
  double weight, rndr, rnd, rnd_xi, rnd_frac, frac1, frac2, ratio;
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
  
  for(i=0; i<parset.n_cloud_per_task; i++)
  {
// generate a direction of the angular momentum of the orbit   
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

    /* Polar coordinates to Cartesian coordinate */
    x = r * cos(phi); 
    y = r * sin(phi);
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
  int i, j, nc;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, Fin, Fout, alpha, k, gam, mu, xi;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, V, Vr, Vph, Vkep, rhoV, theV, linecenter=0.0;
  double cos_Lphi, sin_Lphi, cos_Lthe, sin_Lthe, cos_inc_cmp, sin_inc_cmp;
  double weight, rndr, rnd, rnd_frac, rnd_xi, frac1, frac2, ratio, Rs, g;
  double vx, vy, vz, vxb, vyb, vzb;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model6(const void *pm, int flag_save)
{
  int i, nc;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
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

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  sin_inc_cmp = cos(inc);//sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc);//cos(PI/2.0 - inc);

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

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
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
    /*x = xb * cos(PI/2.0-inc) + zb * sin(PI/2.0-inc);
    y = yb;
    z =-xb * sin(PI/2.0-inc) + zb * cos(PI/2.0-inc); */

    x = xb * cos_inc_cmp + zb * sin_inc_cmp;
    y = yb;
    z =-xb * sin_inc_cmp + zb * cos_inc_cmp;

    dis = r - x;
    weight = 0.5 + k*(x/r);    
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
  int i, j, nc;
  double r, phi, cos_phi, sin_phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb;
  double V, rhoV, theV, Vr, Vph, Vkep, Rs, g;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot, sig_turb;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp, linecenter=0.0;
  double weight, rnd, rnd_xi;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
 * This function caclulate 1d transfer function.
 */
void transfun_1d_cloud_sample_model7(const void *pm, int flag_save)
{
  int i, nc, num_sh;
  double r, phi, dis, Lopn_cos;
  double x, y, z, xb, yb, zb, zb0;
  double inc, F, beta, mu, k, gam, xi, a, s, rin, sig;
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

  a = 1.0/beta/beta;
  s = mu/a;
  rin=mu*F;
  sig=(1.0-F)*s;

  sin_inc_cmp = cos(inc); //sin(PI/2.0 - inc);
  cos_inc_cmp = sin(inc); //cos(PI/2.0 - inc);

  
  /* number of particles in first region */
  
  num_sh = (int)(parset.n_cloud_per_task * model->fsh);

  for(i=0; i<num_sh; i++)
  {
    /* generate a direction of the angular momentum of the orbit */ 
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

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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

  for(i=num_sh; i<parset.n_cloud_per_task; i++)
  {
    /* generate a direction of the angular momentum of the orbit  */ 
    Lphi = 2.0*PI * gsl_rng_uniform(gsl_r);
    Lthe = acos(Lopn_cos_un2 + (Lopn_cos_un1-Lopn_cos_un2) * pow(gsl_rng_uniform(gsl_r), gam));
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

    /* Polar coordinates to Cartesian coordinates */
    x = r * cos(phi); 
    y = r * sin(phi);
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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
  int i, j, nc, num_sh;
  double r, phi, dis, Lopn_cos, cos_phi, sin_phi;
  double x, y, z, xb, yb, zb, zb0, vx, vy, vz, vxb, vyb, vzb, Rs, g, sig_turb;
  double V, rhoV, theV, Vr, Vph, Vkep;
  double inc, F, beta, mu, k, gam, xi, a, s, sig, rin;
  double mbh, fellip, fflow, sigr_circ, sigthe_circ, sigr_rad, sigthe_rad, theta_rot;
  double Lphi, Lthe, sin_Lphi, cos_Lphi, sin_Lthe, cos_Lthe, sin_inc_cmp, cos_inc_cmp,linecenter=0.0;
  double weight, rnd, rnd_xi;
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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;

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

    dis = r - x;
    weight = 0.5 + k*(x/r);
    clouds_tau[i] = dis;
    clouds_weight[i] = weight;


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
  fclose(fp);
}

void restart_action_2d(int iflag)
{
  restart_action_1d(iflag);
}