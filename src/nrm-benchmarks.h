/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <inttypes.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <float.h>

//#include <nrm.h>

#ifndef NRM_BENCHMARKS_H
#define NRM_BENCHMARKS_H 1

#define NRMB_MAX(a, b) (((a) < (b)) ? (a) : (b))
#define NRMB_MIN(a, b) (((a) < (b)) ? (b) : (a))
#define NRMB_ABS(a) (((a) < 0) ? -(a) : (a))

typedef struct timespec nrmb_time_t;

void nrmb_gettime(nrmb_time_t *now);
int64_t nrmb_timediff(const nrmb_time_t *start, const nrmb_time_t *end);

int nrmb_check_double(double ref, double value, int bits);

#endif
