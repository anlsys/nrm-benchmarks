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

#include <nrm.h>
#include <math.h>

#define MK 16
#define NK (1 << MK)
#define NQ 10

/* utils functions */
double randlc (double *, double);
void vranlc (int, double *, double, double []);

static struct nrm_context *context;

static double x[2*NK];
#pragma omp threadprivate(x)
static double q[NQ];

void ep_kernel(double *gc, double *rx, double *ry, double a, double s, double an, size_t nn)
{
	int k_offset = -1;
	double sx = 0.0, sy = 0.0;

#pragma omp parallel copyin(x)
	{
		double t1, t2, t3, t4, x1, x2;
		int k, kk, i, ik, l;
		double qq[NQ];		/* private copy of q[0:NQ-1] */

		for (i = 0; i < NQ; i++) qq[i] = 0.0;

#pragma omp for reduction(+:sx,sy) schedule(static)  
		for (k = 1; k <= (int)nn; k++) {
			kk = k_offset + k;
			t1 = s;
			t2 = an;

			/*      Find starting seed t1 for this kk. */

			for (i = 1; i <= 100; i++) {
				ik = kk / 2;
				if (2 * ik != kk) t3 = randlc(&t1, t2);
				if (ik == 0) break;
				t3 = randlc(&t2, t2);
				kk = ik;
			}

			/*      Compute uniform pseudorandom numbers. */

			vranlc(2*NK, &t1, a, x);

			/*
			   c       Compute Gaussian deviates by acceptance-rejection method and 
			   c       tally counts in concentric square annuli.  This loop is not 
			   c       vectorizable.
			   */

			for ( i = 0; i < NK; i++) {
				x1 = 2.0 * x[2*i] - 1.0;
				x2 = 2.0 * x[2*i+1] - 1.0;
				t1 = pow(x1,2) + pow(x2,2);
				if (t1 <= 1.0) {
					t2 = sqrt(-2.0 * log(t1) / t1);
					t3 = fabs(x1 * t2);				/* Xi */
					t4 = fabs(x2 * t2);				/* Yi */
					l = fmax(t3, t4);
					qq[l] += 1.0;				/* counts */
					sx = sx + t3;				/* sum of Xi */
					sy = sy + t4;				/* sum of Yi */
				}
			}
		}
#pragma omp critical
		{
			for (i = 0; i < NQ; i++) q[i] += qq[i];
		}
	}
	for(size_t i = 0; i < NQ; i++)
		*gc = *gc + q[i];
	*rx = sx;
	*ry = sy;
}





int main(int argc, char **argv)
{
	/* configuration parameters:
	 * - M is the only user-input parameter, it determines the size of the
	 *   problem. All other parameters are derived from it.
	 *   From the NAS NPB documentation:
	 *   - Class: S  W  A  B  C  D  E  F
	 *   - M val: 24 25 28 30 32 36 40 44 
	 */
	size_t m;
	size_t mm, nn;
	double a = 1220703125.0, s = 271828183.0;
	double t, an, gc, rx, ry;
	long int times;

	/* needed for performance measurement */
	int64_t sumtime = 0, mintime = INT64_MAX, maxtime = 0;
	nrmb_time_t start, end;
	int num_threads;

	/* retrieve the size of the problem and initialize the rest of the
	 * configuration from it
	 */
	assert(argc == 3);
	m = strtoull(argv[1], NULL, 0);
	assert(!errno);
	times = strtol(argv[2], NULL, 0);
	assert(!errno);

	mm = m - MK;
	nn = 1 << mm;

	/* ensure that OpenMP is giving us the right number of threads */
#pragma omp parallel
#pragma omp master
	num_threads = omp_get_num_threads();
	int err = 0;
#pragma omp parallel
#pragma omp atomic
	err++;
	assert(num_threads == err);
	err = 0;

	/* initialization: random number generator and private array
	 * dum arrays are there to avoid dead code elimination
	 */
	double dum[3] = { 1.0, 1.0, 1.0 };
	vranlc(0, &(dum[0]), dum[1], &(dum[2]));
	dum[0] = randlc(&(dum[1]), dum[2]);
    
#pragma omp parallel for default(shared)
	for (size_t i = 0; i < 2*NK; i++) x[i] = -1.0e99;

	/* the original NAS EP considers this section part of the benchmark, but
	 * it's single threaded, so we execute it as part of init to make the
	 * whole benchmark even more parallel.
	 */
	vranlc(0, &t, a, x);

	/* Compute AN = A ^ (2 * NK) (mod 2^46). */
	t = a;
	for (size_t i = 0; i <= MK; i++) {
		randlc(&t, t);
	}

	for (size_t i = 0; i < NQ; i++) {
		q[i] = 0.0;
	}
	an = t;
	gc = 0.0;

	/* NRM Context init */
	context = nrm_ctxt_create();
	nrm_init(context, argv[0], 0, 0);

	/* this version of the benchmarks reports one progress each time it goes
	 * through the entire array.
	 */

	nrm_send_progress(context, 1, 1, times, 0, 0, 0);
	for(long int iter = 0; iter < times; iter++)
	{
		int64_t time;
		nrmb_gettime(&start);

		/* the actual benchmark is quite involved,
		 * so we put it in a separate function
		 */
		ep_kernel(&gc, &rx, &ry, a, s, an, nn);

		nrmb_gettime(&end);
		nrm_send_progress(context, 1, 0, 0, 0, 0, 0);

		time = nrmb_timediff(&start, &end);
		sumtime += time;
		mintime = NRMB_MIN(time, mintime);
		maxtime = NRMB_MAX(time, maxtime);
	}

	nrm_fini(context);
	nrm_ctxt_delete(context);

	/* report the configuration and timings */
	fprintf(stdout, "NRM Benchmarks:      %s\n", argv[0]);
	fprintf(stdout, "Version:             %s\n", PACKAGE_VERSION);
	fprintf(stdout, "Description: one progress per iteration, NPB EP benchmark\n");
	fprintf(stdout, "Problem size:        %zu.\n", m);
	fprintf(stdout, "Gaussian Pairs:      %15.0f.\n", gc);
	fprintf(stdout, "Validation values:   %25.15e %25.15e\n", rx, ry);
	fprintf(stdout, "Kernel was executed: %ld times.\n", times);
	fprintf(stdout, "Number of threads:   %d\n", num_threads);
	fprintf(stdout, "Time (s): avg:       %11.6f min:  %11.6f max: %11.6f\n",
		1.0E-09 * sumtime/times, 1.0E-09 * mintime, 1.0E-09 * maxtime);

	/* validate the benchmark: minimum about of bits should be different.
	 * Note that NAS does not give us a validation value for all inputs */
	err = 0;
	switch(m)
	{
		case 24:
			err = err || !nrmb_check_double_prec(1.051299420395306e07, rx, 1e-8);
			err = err || !nrmb_check_double_prec(1.051517131857535e07, ry, 1e-8);
			break;
		case 25:
			err = err || !nrmb_check_double_prec(2.102505525182392e07, rx, 1e-8);
			err = err || !nrmb_check_double_prec(2.103162209578822e07, ry, 1e-8);
			break;
		case 28:
			err = err || !nrmb_check_double_prec(1.682235632304711e08, rx, 1e-8);
			err = err || !nrmb_check_double_prec(1.682195123368299e08, ry, 1e-8);
			break;
		case 30:
			err = err || !nrmb_check_double_prec(6.728927543423024e08, rx, 1e-8);
			err = err || !nrmb_check_double_prec(6.728951822504275e08, ry, 1e-8);
			break;
		case 32:
			err = err || !nrmb_check_double_prec(2.691444083862931e09, rx, 1e-8);
			err = err || !nrmb_check_double_prec(2.691519118724585e09, ry, 1e-8);
			break;
		case 36:
			err = err || !nrmb_check_double_prec(4.306350280812112e10, rx, 1e-8);
			err = err || !nrmb_check_double_prec(4.306347571859157e10, ry, 1e-8);
			break;
		case 40:
			err = err || !nrmb_check_double_prec(6.890169663167274e11, rx, 1e-8);
			err = err || !nrmb_check_double_prec(6.890164670688535e11, ry, 1e-8);
			break;
		case 44:
			err = err || !nrmb_check_double_prec(1.102426773788175e13, rx, 1e-8);
			err = err || !nrmb_check_double_prec(1.102426773787993e13, ry, 1e-8);
			break;
	}
	if(err)
		fprintf(stdout, "VALIDATION FAILED!!!!\n");
	else
		fprintf(stdout, "VALIDATION PASSED!!!!\n");
	return err;
}
