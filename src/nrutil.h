/* BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)ntegrated with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file nrutil.h
 *  \brief Head file for numerical recipes' procedures.
 *
 *  Not used. To be removed.
 */
 
#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void nrerror(char error_text[]);
double *vector(int nl, int nh);
void free_vector(double *v, long nl, long nh);

#endif /* _NR_UTILS_H_ */
