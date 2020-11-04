/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file command_line.h
 *  \brief header file for command_line.c.
 */

#ifndef _BRAINS_COMMANDLINE_H
#define _BRAINS_COMMANDLINE_H

#include <getopt.h>
#include <string.h>

int command_line_options(int argc, char** argv);

#endif