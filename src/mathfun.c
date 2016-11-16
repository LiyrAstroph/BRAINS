/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file mathfun.c
 *  \brief mathematic functions.
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cblas.h>
//#include <clapack.h>
#include <lapacke.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_histogram.h>

#include "allvars.h"
#include "proto.h"

/*!
 * This function calculates matrix multiply C(nxn) = A(nxn) * B(nxn).
 */
void multiply_mat(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply C(nxn) = A^T(nxn) * B(nxn).
 */
void multiply_mat_transposeA(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply C(nxn) = A(nxn) * B^T(nxn).
 */
void multiply_mat_transposeB(double * a, double *b, double *c, int n)
{
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0f
                             , a, n, b, n, 0.0f, c, n);
}

/*!
 * This function calculates matrix multiply Y(n) = A^T(nxn) * X(n).
 */
void multiply_matvec_transposeA(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/*!
 * This function calculates matrix multiply Y(n) = A(nxn) * X(n).
 */
void multiply_matvec(double *a, double *x, int n, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/*!
 * This function calculates  Y(m) = A(m, n) * X(n).
 */
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y)
{
  cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0f, a, n, x, 1, 0.0f, y, 1);
}

/*!
 * This functions calculate A^-1(nxn).
 */
void inverse_mat(double * a, int n, int *info)
{
  int * ipiv;
  ipiv=malloc(n*sizeof(int));

//  dgetrf_(&n, &n, a, &n, ipiv, info);
//  dgetri_(&n, a, &n, ipiv, work, &lwork, info);

  LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, a, n, ipiv);

  free(ipiv);
  return;
}

/*!
 * This function calculate eigenvectors and eigenvalues of matrix.
 */
void eigen_sym_mat(double *a, int n, double *val, int *info)
{
    char jobz='V', uplo='U';

/* store the eigenvectors  in a by rows.
 * store the eigenvalues in val in ascending order.
 */
//    dsyev_(&jobz, &uplo, &n, a, &n, val, work, &lwork, info);

/* store the eigenvectors  in a by columns.
 * store the eigenvalues in val in ascending order.
 */
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, jobz, uplo, n, a, n, val);
    return;
}

/*!
 * This function calculate A(nxn) = X^T(1xn)*X(1xn)
 */
void multiply_vec2mat(double * x, double * a, int n)
{
//  cblas_dsyr(CblasRowMajor, CblasUpper, n, 1.0f, x, 1, a, n);
  int i, j;
  for(i=0; i<n; i++)
    for(j=0; j<=i; j++)
    {
      a[i*n+j] = a[j*n+i] = x[i]*x[j];
    }
}

/*!
 * This function calculates determinant of matrix A. \n
 * There are two versions in the internet, the main difference lies at dealing with the sign of det. \n
 * The version II is verifed to be \b INCORRECT. \n
 * Note that LAPACK is written in C, the indix diffes by 1 with that in C.
 * ********************************
 * Version I:
 * \code{.sh}
 *  det = 1.0; 
 *  for(i=0; i<n; i++)
 *  {
 *   det *= a[i*n+i];
 *   if (ipiv[i] != i+1) 
 *   {
 *     det = -det;
 *   }
 *  }
 * \endcode
 * ********************************
 * Version II:
 * \code{.sh}
 *  det = 1.0;
 *  for(i=0; i<n; i++)
 *  {
 *   det *= a[i*n+i];
 *   if (ipiv[i] != i+1) 
 *   {
 *     ipiv[ipiv[i]-1] = ipiv[i];
 *     det = -det;
 *   }
 *  }
 * \endcode
 * ********************************
 */
double det_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double det;
  ipiv=malloc(n*sizeof(int));
  
  /* LU decomposition */
//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("# Error, Wrong in det_mat!\n");
    exit(-1);
  }

  det = 1.0;
  for(i=0; i<n; i++)
  {
    printf("%d\n", ipiv[i]);
    det *= a[i*n+i];
    if (ipiv[i] != i+1) // note that LAPACK is written in C, the indix diffes by 1 with C.
    {
      det = -det;
    }
  }
  //det=fabs(det);
  free(ipiv);
  return det;
}

/*!
 * This function calculates logarithm determinant of matrix.
 */
double lndet_mat(double *a, int n, int *info)
{
  int * ipiv;
  int i;
  double lndet;
  ipiv=malloc(n*sizeof(int));

  /* LU factorization */
//  dgetrf_(&n, &n, a, &n, ipiv, info);
  *info=LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, a, n, ipiv);
  if(*info!=0)
  {
    printf("Wrong!\n");
    exit(-1);
  }

  lndet = 0.0;
  for(i=0; i<n; i++)
  {
    lndet += log(a[i*n+i]);
  }
  free(ipiv);
  return lndet;
}

/*!
 * This function performs Cholesky decomposition of matrix into upper triangle matrixs
 */
void Chol_decomp_U(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'U';   // decomposite as A = U^T*U, i.e., upper triangle matrix.
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
//    exit(-1);
    return;
  }
  else if (*info>0)
  {
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
//    exit(-1);
    return;
  }

  // only the upper triangle is referenced by dpotrf output, 
  // so the strictly lower triangle are must set to zero
  for(i=0;i<n;i++)
    for(j=0;j<i;j++)
      a[i*n+j] = 0.0;
  return;
}

/*!
 * This function performs Cholesky decomposition of matrix into lower triangle matrixs
 */
void Chol_decomp_L(double *a, int n, int *info)
{
  int i,j;
  char uplo = 'L';   // decomposite as A = L*L^T, i.e., lower triangle matrix.
//  dpotrf_(&uplo, &n, a, &n, info);
  *info=LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, n);
  if(*info<0)
  {
    fprintf(stderr, "The %d-th argument had an illegal value!\n", *info);
//    exit(-1);
    return;
  }
  else if (*info>0)
  {
    fprintf(stderr, "The leading minor of order %d is not positive definite, and the factorization could not be completed.\n", *info);
//    exit(-1);
    return;
  }

  // only the lower triangle is referenced by dpotrf output, 
  // so the strictly upper triangle are must set to zero
  for(i=0;i<n;i++)
    for(j=i+1;j<n;j++)
      a[i*n+j] = 0.0;
  return;
}

/*!
 * This function display matrix on the screen.
 */
void display_mat(double *a, int m, int n)
{
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("%e\t", a[i*n+j]);
        }
        printf("\n");
    }
}

/*!
 * This function allocates memory for matrix.
 */
double ** matrix_malloc(int n1, int n2)
{
  double ** mat;
  int i;

  if(!(mat = malloc(n1*sizeof(double*))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  for(i=0; i<n1; i++)
  {
    if(!(mat[i] = malloc(n2*sizeof(double))))
    {
      fprintf(stderr, "Unable to allocate the matrix!\n");
      exit(-1);
    }
  }
  return mat;
}

/*!
 * This function allocates memory for array.
 */
double * array_malloc(int n)
{
  double *array;

  if(!(array = malloc(n*sizeof(double))))
  {
    fprintf(stderr, "Unable to allocate the matrix!\n");
    exit(-1);
  }

  return array;
}

/*!
 * This function is to test the functions defined in mathfun.c.
 */
void test_mathfun()
{
  double *A, det;
  int n = 3, info;
  A = malloc(n*n*sizeof(double));

  A[0*n+0] = 1.0;
  A[0*n+1] = 10.3;
  A[0*n+2] = 6.3;

  A[1*n+0] = -8.3;
  A[1*n+1] = -8.0;
  A[1*n+2] = 5.3;

  A[2*n+0] = 0.3;
  A[2*n+1] = -9.3;
  A[2*n+2] = -3.0;
  
  det = det_mat(A, n, &info);
  printf("%f\n", det);
}