/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

#ifndef _VERSION_H
#define _VERSION_H

/* BRAINS version number */
#define BRAINS_MAJOR_VERSION 0
#define BRAINS_MINOR_VERSION 1
#define BRAINS_PATCH_VERSION 5

void print_version();
void fprint_version(FILE *fp);
#endif