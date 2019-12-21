/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 * \file user_transfun.c
 * \brief allow users to define their own transfer function.
 * 
 * User should provide: 1) number of parameters in 1D and 2D; 
 *                      2) function for setting parameter range
 *                      3) function for calculating 1D transfer function
 *                      4) function for calculating 2D transfer function
 *                      5) function for setting the parameter values used in simulation (if flag_dim<0)
 *
 * The following is an example. The users can make modification with their own BLR models.
 *
 */

#include "brains.h"

const int num_params_MyTransfun1d = 2; /* number of parameters for 1D transfer function */
const int num_params_MyTransfun2d = 4; /* number of parameters for 2D transfer function */

/*!
 * this function set the parameter prior range.
 */
void set_par_range_mytransfun()
{
  int i;

  i=0;
  blr_range_model[i][0] = 0.0;
  blr_range_model[i++][1] = rcloud_max_set*0.5;

  blr_range_model[i][0] = log(0.1);
  blr_range_model[i++][1] = log(100.0);

  return;
}

/*!
 * this function calculates 1D transfer function.
 */
void transfun_1d_cal_mytransfun(const void *pm, int flag_save)
{
  int i;
  double cen, sig, Anorm;
  double tau_min, tau_max;
  double *model = (double *)pm;
  
  cen = model[0];
  sig = exp(model[1]);
  
  tau_min = 0.0;
  tau_max = cen + sig * 10.0;
  dTransTau = (tau_max - tau_min)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = tau_min + dTransTau * i;
    Trans1D[i] = exp( -0.5 * pow((TransTau[i] - cen)/sig, 2.0) );
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
 * this function calculates 2D transfer function.
 */
void transfun_2d_cal_mytransfun(const void *pm, double *transv, double *trans2d, int n_vel, int flag_save)
{
  int i, j;
  double Anorm, dV;
  double cen, sig, cenV, sigV;
  double tau_min, tau_max;
  double *model = (double *)pm;

  cen = model[0];
  sig = exp(model[1]);
  cenV = model[2];
  sigV = exp(model[3]);
  
  tau_min = 0.0;
  tau_max = cen + sig * 10.0;
  dTransTau = (tau_max - tau_min)/(parset.n_tau - 1);
  for(i=0; i<parset.n_tau; i++)
  {
    TransTau[i] = tau_min + dTransTau * i;
  }

  dV =(transv[1] - transv[0]); /* velocity grid width */
  for(i=0; i<parset.n_tau; i++)
    for(j=0;j<n_vel;j++)
      trans2d[i*n_vel+j]= exp( -0.5 * pow((TransTau[i] - cen)/sig, 2.0) ) 
                        * exp( -0.5 * pow((transv[j] - cenV)/sigV, 2.0) ) ; 

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

/*!
 * this function set parameter value used in simulations.
 */
void set_par_value_mytransfun_sim(double *pm)
{
  int i;

	i=0;
  pm[i++] = 4.0;
  pm[i++] = log(1.0);
  pm[i++] = 0.0;
  pm[i++] = log(1000.0/VelUnit);
  return;
}