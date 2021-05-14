/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

#include "config.h"

#include "nrm-benchmarks.h"

void nrmb_gettime(nrmb_time_t *now)
{
	clock_gettime(CLOCK_REALTIME, now);
}

int64_t nrmb_timediff(const nrmb_time_t *start, const nrmb_time_t *end)
{
	int64_t timediff = (end->tv_nsec - start->tv_nsec) +
		1000000000 * (end->tv_sec - start->tv_sec);
	return timediff;
}

int nrmb_check_double(double ref, double value, int bits)
{
	double diff = NRMB_ABS(ref - value);
	return diff <= NRMB_MAX(NRMB_ABS(ref), NRMB_ABS(value)) * DBL_EPSILON *
		((1 << bits) - 1);
}

int nrmb_check_double_prec(double ref, double value, double eps)
{
	double diff = NRMB_ABS(ref - value);
	return diff <= NRMB_MAX(NRMB_ABS(ref), NRMB_ABS(value)) * eps;
}
