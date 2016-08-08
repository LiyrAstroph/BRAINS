#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)


#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
float *vector(int nl, int nh);
void free_vector(float *v, long nl, long nh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
void free_vector();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */
