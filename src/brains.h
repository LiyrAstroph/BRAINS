/*
 * BRAINS
 * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
 * Yan-Rong Li, liyanrong@ihep.ac.cn
 * Thu, Aug 4, 2016
 */

/*!
 *  \file brains.h
 *  \brief include all the header files.
 */

#ifndef _BRAINS_H
#define _BRAINS_H

#include "allvars.h"
#include "proto.h"
#include "blr_models.h"
#include "version.h"
#include "dnest_con.h"
#include "dnest_line1d.h"
#include "dnest_line2d.h"
#include "command_line.h"
#include "user_blr_model.h"
#include "user_transfun.h"

#ifdef SA
#include "dnest_sa.h"
#include "dnest_sa1d.h"
#include "dnest_sa2d.h"
#endif

/* header files for CDNest */
#include "dnestvars.h"

#endif