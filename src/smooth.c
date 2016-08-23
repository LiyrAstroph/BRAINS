/*
 * BRAINS
 * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

#define NRANSI
 
#include "nrutil.h"

void convlv(float data[], unsigned long n, float respns[], unsigned long m,
  int isign, float ans[]);

/* Gaussian smooth */
void line_gaussian_smooth_2D_FFT(double *transv, double *fl2d, int nl, int nv)
{
  int i, j, nd, nM, isign;
  double sigV;
  float *data, *respns, *resp, *ans, tot, dV;

  sigV = parset.InstRes / VelUnit;
  nd = log(nv)/log(2.0) + 1.0;
  nd = 1<<nd;
  nM = (nd/2) *2 - 1;

  data = vector(1, nd);
  respns = vector(1, nd);
  resp = vector(1, nd);
  ans = vector(1, nd*2);
  
  dV = TransV[1] - TransV[0];
  tot = 0.0;
  for(i=1; i<=nM/2; i++)
  {
    respns[i] = exp(-((i-1)*dV)*((i-1)*dV)/sigV/sigV/2.0);
    tot += respns[i];
  }

  for(i=nM; i>=nM/2+2; i--)
  {
    respns[i] = exp(-((nM-i+1)*dV)*((nM-i+1)*dV)/sigV/sigV/2.0);
    tot += respns[i];
  }

  i = nM/2+1;
  respns[nM/2+1] = exp(-((i-1)*dV)*((i-1)*dV)/sigV/sigV/2.0);
  tot += respns[nM/2+1];
  
  for(i=1; i<=nM; i++)
  {
    respns[i] /= (tot);
  }

  for(j=0; j<nl; j++)
  {
    memcpy(resp+1, respns+1, nd*sizeof(float));
    for(i=0; i<nv; i++)
    {
      data[i+1] = fl2d[j*nv + i];
    }
    for(i=nv; i<nd; i++)
      data[i+1] = 0.0;

    isign = 1;
    convlv(data,nd,resp,nM,isign,ans);

    for(i=0; i<nv; i++)
    {
      fl2d[j*nv + i] = ans[i+1];      
//      printf("%f %f %f %f\n", TransV[i]*VelUnit, ans[i+1], data[i+1], respns[i+1]);
    }
  }

  free_vector(ans, 1, 2*nd);
  free_vector(respns, 1, nd);
  free_vector(data, 1, nd);
  free_vector(resp, 1, nd);

}

void convlv(float data[], unsigned long n, float respns[], unsigned long m,
	int isign, float ans[])
{
	void realft(float data[], unsigned long n, int isign);
	void twofft(float data1[], float data2[], float fft1[], float fft2[],
		unsigned long n);
	unsigned long i,no2;
	float dum,mag2,*fft;

	fft=vector(1,n<<1);
	for (i=1;i<=(m-1)/2;i++)
		respns[n+1-i]=respns[m+1-i];
	for (i=(m+3)/2;i<=n-(m-1)/2;i++)
		respns[i]=0.0;
	twofft(data,respns,fft,ans,n);
	no2=n>>1;
	for (i=2;i<=n+2;i+=2) {
		if (isign == 1) {
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
			ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
		} else if (isign == -1) {
			if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
				nrerror("Deconvolving at response zero in convlv");
			ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;
			ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
		} else nrerror("No meaning for isign in convlv");
	}
	ans[2]=ans[n+1];
	realft(ans,n,-1);
	free_vector(fft,1,n<<1);
}

#undef NRANSI

void realft(float data[], unsigned long n, int isign)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long i,i1,i2,i3,i4,np3;
	float c1=0.5,c2,h1r,h1i,h2r,h2i;
	double wr,wi,wpr,wpi,wtemp,theta;

	theta=3.141592653589793/(double) (n>>1);
	if (isign == 1) {
		c2 = -0.5;
		four1(data,n>>1,1);
	} else {
		c2=0.5;
		theta = -theta;
	}
	wtemp=sin(0.5*theta);
	wpr = -2.0*wtemp*wtemp;
	wpi=sin(theta);
	wr=1.0+wpr;
	wi=wpi;
	np3=n+3;
	for (i=2;i<=(n>>2);i++) {
		i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
		h1r=c1*(data[i1]+data[i3]);
		h1i=c1*(data[i2]-data[i4]);
		h2r = -c2*(data[i2]+data[i4]);
		h2i=c2*(data[i1]-data[i3]);
		data[i1]=h1r+wr*h2r-wi*h2i;
		data[i2]=h1i+wr*h2i+wi*h2r;
		data[i3]=h1r-wr*h2r+wi*h2i;
		data[i4] = -h1i+wr*h2i+wi*h2r;
		wr=(wtemp=wr)*wpr-wi*wpi+wr;
		wi=wi*wpr+wtemp*wpi+wi;
	}
	if (isign == 1) {
		data[1] = (h1r=data[1])+data[2];
		data[2] = h1r-data[2];
	} else {
		data[1]=c1*((h1r=data[1])+data[2]);
		data[2]=c1*(h1r-data[2]);
		four1(data,n>>1,-1);
	}
}

void twofft(float data1[], float data2[], float fft1[], float fft2[],
	unsigned long n)
{
	void four1(float data[], unsigned long nn, int isign);
	unsigned long nn3,nn2,jj,j;
	float rep,rem,aip,aim;

	nn3=1+(nn2=2+n+n);
	for (j=1,jj=2;j<=n;j++,jj+=2) {
		fft1[jj-1]=data1[j];
		fft1[jj]=data2[j];
	}
	four1(fft1,n,1);
	fft2[1]=fft1[2];
	fft1[2]=fft2[2]=0.0;
	for (j=3;j<=n+1;j+=2) {
		rep=0.5*(fft1[j]+fft1[nn2-j]);
		rem=0.5*(fft1[j]-fft1[nn2-j]);
		aip=0.5*(fft1[j+1]+fft1[nn3-j]);
		aim=0.5*(fft1[j+1]-fft1[nn3-j]);
		fft1[j]=rep;
		fft1[j+1]=aim;
		fft1[nn2-j]=rep;
		fft1[nn3-j] = -aim;
		fft2[j]=aip;
		fft2[j+1] = -rem;
		fft2[nn2-j]=aip;
		fft2[nn3-j]=rem;
	}
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void four1(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
#undef SWAP