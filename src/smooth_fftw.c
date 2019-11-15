/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file smooth.c
 *  \brief smoothing functions using FFT provided by FFTW library.
 *
 *  References: Press et al., Numerical Recipes, Chapter 13.1.
 *              FFTW website: http://www.fftw.org/
 *
 *  To facilitate computation, we use the factor:
 *  the Fourier transform a Gaussian is still a Gaussian, which can be expressed analytically.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>
#include <fftw3.h>

#include "brains.h"


int nd_fft, npad;

fftw_complex *data_fft, *resp_fft, *conv_fft;
fftw_plan pdata, presp, pback;
double *real_data, *real_resp, *real_conv;


/*!
 * This function initates workspace for FFT.
 */
void smooth_init(int nv, const double *transv)
{
  npad = fmin(nv * 0.2, 100);
  npad = (npad/2) * 2;
  nd_fft = nv+npad;

  data_fft = (fftw_complex *) fftw_malloc((nd_fft/2+1) * sizeof(fftw_complex));
  resp_fft = (fftw_complex *) fftw_malloc((nd_fft/2+1) * sizeof(fftw_complex));
  conv_fft = (fftw_complex *) fftw_malloc((nd_fft/2+1) * sizeof(fftw_complex));

  real_data = (double *)fftw_malloc(nd_fft * sizeof(double));
  //real_resp = (double *)fftw_malloc(nd_fft * sizeof(double));
  real_conv = (double *)fftw_malloc(nd_fft * sizeof(double));

  pdata = fftw_plan_dft_r2c_1d(nd_fft, real_data, data_fft, FFTW_MEASURE);
  //presp = fftw_plan_dft_r2c_1d(nd_fft, real_resp, resp_fft, FFTW_MEASURE);
  pback = fftw_plan_dft_c2r_1d(nd_fft, conv_fft, real_conv, FFTW_MEASURE);

}

/*!
 * This function finalizes FFT.
 */
void smooth_end()
{
  fftw_destroy_plan(pdata);
  //fftw_destroy_plan(presp);
  fftw_destroy_plan(pback);

  fftw_free(data_fft);
  fftw_free(resp_fft);
  fftw_free(conv_fft);
  
  fftw_free(real_data);
  //fftw_free(real_resp);
  fftw_free(real_conv);
}

/*!
 * This function performs FFT-smoothing to 2d line.
 */
void line_gaussian_smooth_2D_FFT(const double *transv, double *fl2d, int nl, int nv, const void *pm)
{
  int i, j;
  // initialize response and its fft.
  double sigV, sigV_instrinsic, dV;
  double *pmodel = (double *)pm;

  dV = transv[1] - transv[0];

  if(parset.flag_narrowline < 2)
  {
    sigV_instrinsic = parset.width_narrowline;
  }
  else
  {
    sigV_instrinsic = parset.width_narrowline + pmodel[num_params_blr-num_params_res-num_params_linecenter-1-2] * parset.width_narrowline_err;
  }
  
  if(parset.InstRes >= 0.0)
  {
    sigV = parset.InstRes + pmodel[num_params_blr-num_params_res-num_params_linecenter-1]*parset.InstRes_err;
    sigV = fmax(0.0, sigV*sigV - sigV_instrinsic*sigV_instrinsic);
    sigV = sqrt(sigV);

    /* setup response 
       note the factor nd_fft. 
       FFTW does not include the factor 1/nd_fft when transforming backforward.
       So include here. Thereby, need not to mannually multiply 1/nd_fft in the final output.
     */

    for(i=0; i<nd_fft/2+1; i++)
    {
      resp_fft[i][0] = exp(-2.0 * PI*PI * sigV/dV*sigV/dV * i*i*1.0/nd_fft/nd_fft)/nd_fft;
      resp_fft[i][1] = 0.0;
    }
    
    for(j=0; j<nl; j++)
    {
      memcpy(real_data+npad/2, &fl2d[j*nv], nv*sizeof(double));
      for(i=0; i<npad/2; i++)
        real_data[i] = real_data[nd_fft-1-i] = 0.0;

      /* FFT of line */
      fftw_execute(pdata);
    
      /* complex multiply and inverse FFT 
       * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
       * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
       */
      for(i=0; i<nd_fft/2 + 1; i++)
      {
        conv_fft[i][0] = data_fft[i][0]*resp_fft[i][0] - data_fft[i][1]*resp_fft[i][1];
        conv_fft[i][1] = data_fft[i][0]*resp_fft[i][1] + data_fft[i][1]*resp_fft[i][0];
      }
      fftw_execute(pback);

      //for(i=0; i<nv; i++)
      //{
      //fl2d[j*nv + i] = real_conv[i] * dV / nd_fft;
      //}
      memcpy(&fl2d[j*nv], real_conv+npad/2, nv*sizeof(double));
    }

  }
  else
  {
    for(j=0; j<nl; j++)
    {
      sigV = instres_epoch[j] + pmodel[num_params_blr-num_params_res-num_params_linecenter-1 + j]*instres_err_epoch[j];
      sigV = fmax(0.0, sigV*sigV - sigV_instrinsic*sigV_instrinsic);
      sigV = sqrt(sigV);

      /* setup response */
      for(i=0; i<nd_fft/2+1; i++)
      {
        resp_fft[i][0] = exp(-2.0 * PI*PI * sigV/dV*sigV/dV * i*i*1.0/nd_fft/nd_fft)/nd_fft;
        resp_fft[i][1] = 0.0;
      }
    
      memcpy(real_data+npad/2, &fl2d[j*nv], nv*sizeof(double));
      for(i=0; i<npad/2; i++)
        real_data[i] = real_data[nd_fft-1-i] = 0.0;

      /* FFT of line */
      fftw_execute(pdata);
    
      /* complex multiply and inverse FFT 
       * note that for FFT of real data, FFTW outputs n/2+1 complex numbers.
       * similarly, for complex to real transform, FFTW needs input of n/2+1 complex numbers.
       */
      for(i=0; i<nd_fft/2 + 1; i++)
      {
        conv_fft[i][0] = data_fft[i][0]*resp_fft[i][0] - data_fft[i][1]*resp_fft[i][1];
        conv_fft[i][1] = data_fft[i][0]*resp_fft[i][1] + data_fft[i][1]*resp_fft[i][0];
      }
      fftw_execute(pback);

      //for(i=0; i<nv; i++)
      //{
        //fl2d[j*nv + i] = real_conv[i] * dV / nd_fft;
      //}
      memcpy(&fl2d[j*nv], real_conv+npad/2, nv*sizeof(double));
    }
  }
  return;
}


/*!
 * convolution test.
 */
void smooth_test()
{
  int i, j;
  double sig_a, sig_b;
  double dV, fb, fa;

  const int N=100;
  fftw_complex *out_a, *out_b, *out;
  double *in, *in_a, *re, *in_b;
  fftw_plan pa, pb, pout;
  
  in = (double*) fftw_malloc(sizeof(double) * N);
  in_a = (double*) fftw_malloc(sizeof(double) * N);
  in_b = (double*) fftw_malloc(sizeof(double) * N);
  re = (double*) fftw_malloc(sizeof(double) * N);
  
  out_a = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
  out_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));
  
  pa = fftw_plan_dft_r2c_1d(N, in_a, out_a, FFTW_MEASURE);
  pb = fftw_plan_dft_r2c_1d(N, in_b, out_b, FFTW_MEASURE); 
  pout = fftw_plan_dft_c2r_1d(N, out, re, FFTW_MEASURE);
  
  dV = 2.0;
  
  fb = 5.0;
  fa = 10.0;
  
  for(j=0; j<1; j++)
  {
  
    sig_a = 10.0*(j+1);
    sig_b = 5.0*(j+1);
  
    for(i=0; i<N; i++)
    {
      in_a[i] = 1.0/sqrt(2.0*PI)/sig_a * exp(-0.5 * (i-N/2.0)*(i-N/2.0)*dV*dV/sig_a/sig_a);
      in_a[i] = in_a[i]/fa;
    }
    memcpy(in, in_a, N*sizeof(double));
  
    fftw_execute(pa);
  
    /*tot = 0.0;
    for(i=0; i<N/2; i++)
    {
      in_b[i] = 1.0/sqrt(2.0*PI)/sig_b * exp(-0.5 * (i*dV)*(i*dV)/sig_b/sig_b);
      tot += in_b[i];
    }
    for(i=N-1; i>N/2; i--)
    {
      in_b[i] = 1.0/sqrt(2.0*PI)/sig_b * exp(-0.5 * (i-N)*(i-N)*dV*dV/sig_b/sig_b);
      tot += in_b[i];
    }
  
    for(i=0; i<N; i++)
    {
      in_b[i] = in_b[i]/ (tot) / N / fb;
    }
 
    fftw_execute(pb);*/

    for(i=0; i<N/2+1; i++)
    {
      out_b[i][0] = exp(-2.0 * PI*PI * sig_b/dV*sig_b/dV * i*i*1.0/N/N)/N/fb;
      out_b[i][1] = 0.0; //-sin(2.0 * PI*PI * sig_b/dV*sig_b/dV * i*i*1.0/N/N)/N/fb;
    }
    
    /*for(i=0; i<N/2+1; i++)
    {
      printf("%f %f %f\n", i*1.0/N, out_b[i][0], out_b[i][1]);
    }*/

  
    for(i=0; i<N/2+1; i++)
    {
      out[i][0] = out_a[i][0]*out_b[i][0] - out_a[i][1]*out_b[i][1];
      out[i][1] = out_a[i][0]*out_b[i][1] + out_a[i][1]*out_b[i][0];
    }
  
    fftw_execute(pout);
  
    for(i=0; i<N; i++)
    {
      printf("%d %f %f %f\n", i, in[i],  re[i] * fb, 
        1.0/sqrt(2.0*PI)/sqrt(sig_a*sig_a+sig_b*sig_b) * exp(-0.5 * (i-N/2.0)*(i-N/2.0)*dV*dV/(sig_a*sig_a+sig_b*sig_b))/fa);
    }
  }
  
  fftw_destroy_plan(pa);
  fftw_destroy_plan(pb);
  fftw_destroy_plan(pout);
  fftw_free(out_a);
  fftw_free(out_b);
  fftw_free(out);
  fftw_free(in_a);
  fftw_free(in_b);
  fftw_free(in);
  fftw_free(re);
}
