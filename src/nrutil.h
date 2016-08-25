#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

void nrerror(char error_text[]);
double *vector(int nl, int nh);
void free_vector(double *v, long nl, long nh);

#endif /* _NR_UTILS_H_ */
