/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file smooth.c
 *  \brief smoothing functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_halfcomplex.h>

#include "allvars.h"
#include "proto.h"

int nd_fft;

double *resp, *data_fft;
double *resp_cmp, *data_fft_cmp, *data_fft_inverse;

gsl_fft_real_wavetable * real_data, *real_resp;
gsl_fft_complex_wavetable * cmp_data;
gsl_fft_real_workspace * work_data, *work_resp;
gsl_fft_complex_workspace * work_cmp;

void smooth_init(int nv)
{
  //nd_fft = parset.n_vel_recon>=n_vel_data? parset.n_vel_recon : n_vel_data;
  nd_fft = nv;
  //printf("%d\n", nd_fft);
  real_data = gsl_fft_real_wavetable_alloc (nd_fft);
  real_resp = gsl_fft_real_wavetable_alloc (nd_fft);
  cmp_data = gsl_fft_complex_wavetable_alloc (nd_fft);

  work_data = gsl_fft_real_workspace_alloc (nd_fft);
  work_resp = gsl_fft_real_workspace_alloc (nd_fft);
  work_cmp = gsl_fft_complex_workspace_alloc (nd_fft);

  resp = malloc(nd_fft * sizeof(double));
  data_fft = malloc(nd_fft * sizeof(double));

  resp_cmp = malloc(nd_fft * 2 * sizeof(double));
  data_fft_cmp = malloc(nd_fft *2 * sizeof(double));
  data_fft_inverse = malloc(nd_fft *2 * sizeof(double));
}

void smooth_end()
{

  gsl_fft_real_wavetable_free(real_data);
  gsl_fft_real_wavetable_free(real_resp);
  gsl_fft_complex_wavetable_free(cmp_data);
  
  gsl_fft_real_workspace_free(work_data);
  gsl_fft_real_workspace_free(work_resp);
  gsl_fft_complex_workspace_free(work_cmp);
  
  free(resp);
  free(data_fft);

  free(resp_cmp);
  free(data_fft_cmp);
  free(data_fft_inverse);

}
void line_gaussian_smooth_2D_FFT(const double *transv, double *fl2d, int nl, int nv)
{
  int i, j;
  double sigV, dV, tot;

  sigV = parset.InstRes / VelUnit;
  dV = transv[1] - transv[0];

  tot = 0.0;
  for (i = 0; i<nd_fft/2; i++)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/sigV * exp(-0.5*(i*dV)*(i*dV)/sigV/sigV);
    tot += resp[i];
  }
  for (i = nd_fft-1; i>= nd_fft/2; i--)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/sigV * exp(-0.5*((i-nd_fft)*dV)*((i-nd_fft)*dV)/sigV/sigV);
    tot += resp[i];
  }  
  
  for(i=0; i<nd_fft; i++)
  {
    resp[i] /= (tot * dV);
  }
  gsl_fft_real_transform(resp, 1, nd_fft, real_resp, work_resp);
  gsl_fft_halfcomplex_unpack(resp, resp_cmp, 1, nd_fft);

  for(j=0; j<nl; j++)
  {
    memcpy(data_fft, &fl2d[j*nv], nv*sizeof(double));

    gsl_fft_real_transform(data_fft, 1, nd_fft, real_data, work_data);
    gsl_fft_halfcomplex_unpack(data_fft, data_fft_cmp, 1, nd_fft);
    
    for(i=0; i<nd_fft; i++)
    {
      data_fft_inverse[i*2] = data_fft_cmp[i*2]*resp_cmp[i*2] - data_fft_cmp[i*2+1]*resp_cmp[i*2+1];
      data_fft_inverse[i*2+1] = data_fft_cmp[i*2]*resp_cmp[i*2+1] + data_fft_cmp[i*2+1]*resp_cmp[i*2];
    }
    gsl_fft_complex_inverse(data_fft_inverse, 1, nd_fft, cmp_data, work_cmp);

    for(i=0; i<nv; i++)
    {
      // take into account the velocity grid width
      fl2d[j*nv + i] = data_fft_inverse[i*2] * dV;
    }
  }
  return;
}

// test smooth 
void line_gaussian_smooth_2D_FFT_test(const double *transv, double *fl2d, int nl, int nv)
{
  int i, j;
  double sigV, dV, tot;

  sigV = parset.InstRes / VelUnit;
  dV = transv[1] - transv[0];

  tot = 0.0;
  for (i = 0; i<nd_fft/2; i++)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/sigV * exp(-0.5*(i*dV)*(i*dV)/sigV/sigV);
    tot += resp[i];
  }
  for (i = nd_fft-1; i>= nd_fft/2; i--)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/sigV * exp(-0.5*((i-nd_fft)*dV)*((i-nd_fft)*dV)/sigV/sigV);
    tot += resp[i];
  }  
  
  for(i=0; i<nd_fft; i++)
  {
    resp[i] /= (tot * dV);
  }
  gsl_fft_real_transform(resp, 1, nd_fft, real_resp, work_resp);
  gsl_fft_halfcomplex_unpack(resp, resp_cmp, 1, nd_fft);

  for(j=0; j<nl; j++)
  {
    //memcpy(data_fft, &fl2d[j*nv], nv*sizeof(double));
    for(i=0; i<nv; i++)
    {
      data_fft[i] =   exp( - 0.5 * pow(transv[i]/(1500.0/VelUnit), 2.0));
    }

    gsl_fft_real_transform(data_fft, 1, nd_fft, real_data, work_data);
    gsl_fft_halfcomplex_unpack(data_fft, data_fft_cmp, 1, nd_fft);
    
    for(i=0; i<nd_fft; i++)
    {
      data_fft_inverse[i*2] = data_fft_cmp[i*2]*resp_cmp[i*2] - data_fft_cmp[i*2+1]*resp_cmp[i*2+1];
      data_fft_inverse[i*2+1] = data_fft_cmp[i*2]*resp_cmp[i*2+1] + data_fft_cmp[i*2+1]*resp_cmp[i*2];
    }
    gsl_fft_complex_inverse(data_fft_inverse, 1, nd_fft, cmp_data, work_cmp);

    for(i=0; i<nv; i++)
    {
      // take into account the velocity grid width
      fl2d[j*nv + i] = data_fft_inverse[i*2] * dV;
    }
  }
  return;
}

void smooth_test()
{
  int i, n = 205;
  double data[n], resp[n];
  double data_cmp[n*2], resp_cmp[n*2];
  double ans[n*2];
  
  gsl_fft_real_wavetable * real, *real_resp;
  gsl_fft_complex_wavetable * hc;
  gsl_fft_real_workspace * work, *work_resp;
  gsl_fft_complex_workspace * work_cmp;
  
  for (i = 0; i < n; i++)
  {
    data[i] = 1.0/sqrt(2.0*M_PI)/6.0 * exp(-0.5*(i-50.0)*(i-50.0)/6.0/6.0);
  }
  
  for (i = 0; i < n/2; i++)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/10.0 * exp(-0.5*(i-0.0)*(i-0.0)/10.0/10.0);
  }
  for (i = n-1; i >= n/2; i--)
  {
    resp[i] = 1.0/sqrt(2.0*M_PI)/10.0 * exp(-0.5*(i-n)*(i-n)/10.0/10.0);
  }  

  for (i = 0; i < n; i++)
  {
    printf ("%d: %e\n", i, data[i]);
  }
  printf ("\n");

  work = gsl_fft_real_workspace_alloc (n);
  real = gsl_fft_real_wavetable_alloc (n);
  
  work_resp = gsl_fft_real_workspace_alloc (n);
  real_resp = gsl_fft_real_wavetable_alloc (n);

  work_cmp = gsl_fft_complex_workspace_alloc(n);
  hc = gsl_fft_complex_wavetable_alloc (n);

  gsl_fft_real_transform (data, 1, n, real, work);
  gsl_fft_real_transform (resp, 1, n, real_resp, work_resp);
    
  gsl_fft_halfcomplex_unpack(data, data_cmp, 1, n);
  gsl_fft_halfcomplex_unpack(resp, resp_cmp, 1, n);
  
  for(i=0; i<n; i++)
  {
    ans[i*2] = data_cmp[i*2]*resp_cmp[i*2] - data_cmp[i*2+1]*resp_cmp[i*2+1];
    ans[i*2+1] = data_cmp[i*2]*resp_cmp[i*2+1] + data_cmp[i*2+1]*resp_cmp[i*2];
  }

  gsl_fft_complex_inverse(ans, 1, n, hc, work_cmp);
  

  for (i = 0; i < n; i++)
  {
    printf ("%d: %e %e\n", i,  ans[i*2], ans[i*2+1]);
  }
  
  gsl_fft_real_wavetable_free (real);
  gsl_fft_real_wavetable_free (real_resp);
  gsl_fft_real_workspace_free (work);
  gsl_fft_real_workspace_free (work_resp);
  gsl_fft_complex_wavetable_free (hc);
  gsl_fft_complex_workspace_free(work_cmp);
  return;
}