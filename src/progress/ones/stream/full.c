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

static double *a, *b, *c;
static struct nrm_context *context;

int main(int argc, char **argv)
{
	/* configuration parameters:
	 * - array size in number of double elements
	 * - number of times to run through the benchmark
	 */
	size_t array_size;
	long int times;
	double scalar = 3.0;

	/* needed for performance measurement */
	int64_t sumtime[4] = {0,0,0,0};
	int64_t mintime[4] = {INT64_MAX, INT64_MAX, INT64_MAX, INT64_MAX};
	int64_t maxtime[4] = {0, 0, 0, 0};
	const char *names[4] = {"Copy", "Scale", "Add", "Triad"};
	size_t bytes[4] = {2, 2, 3, 3};
	nrmb_time_t progress_start, progress_end;
	int64_t progress_time;
	nrmb_time_t start, end;
	size_t memory_size;
	int num_threads;

	/* retrieve the size of the allocation and the number of time
	 * to loop through the kernel.
	 */
	assert(argc == 3);
	errno = 0;
	array_size = strtoull(argv[1], NULL, 0);
	assert(!errno);
	errno = 0;
	times = strtol(argv[2], NULL, 0);
	assert(!errno);

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

	/* allocate the arrays and initialize them. Note that we expect the
	 * first-touch policy of Linux to result in the arrays being properly
	 * balanced between threads/numa-nodes
	 */
	memory_size = array_size * sizeof(double);
	a = malloc(memory_size);
	b = malloc(memory_size);
	c = malloc(memory_size);

#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
	{
		a[i] = 1.0;
		b[i] = 2.0;
		c[i] = 0.0;
	}

	/* NRM Context init */
	context = nrm_ctxt_create();
	nrm_init(context, argv[0], 0, 0);
	nrmb_gettime(&progress_start);

	/* one run of the benchmark for free, warms up the memory */
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		c[i] = a[i];
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		b[i] = scalar*c[i];
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		c[i] = a[i] + b[i];
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		a[i] = b[i] + scalar*c[i];

	/* this version of the benchmarks reports one progress each time it goes
	 * through the entire array.
	 */
	nrm_send_progress(context, 1);

	for(long int iter = 0; iter < times; iter++)
	{
		int64_t time;

#define TSTART(k) nrmb_gettime(&start)
#define TEND(i) do { \
		nrmb_gettime(&end); \
		time = nrmb_timediff(&start, &end); \
		sumtime[i] += time; \
		mintime[i] = NRMB_MIN(time, mintime[i]); \
		maxtime[i] = NRMB_MAX(time, maxtime[i]); \
	} while(0)

		/* the actual benchmarks */
		TSTART(0);
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		c[i] = a[i];
		TEND(0);
		TSTART(1);
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		b[i] = scalar*c[i];
		TEND(1);
		TSTART(2);
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		c[i] = a[i] + b[i];
		TEND(2);
		TSTART(3);
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		a[i] = b[i] + scalar*c[i];
		TEND(3);

		nrm_send_progress(context, 1);
	}

	nrmb_gettime(&progress_end);
	nrm_fini(context);
	nrm_ctxt_delete(context);
	progress_time = nrmb_timediff(&progress_start, &progress_end);
	/* compute stats */

	/* report the configuration and timings */
	fprintf(stdout, "NRM Benchmarks:      %s\n", argv[0]);
	fprintf(stdout, "Version:             %s\n", PACKAGE_VERSION);
	fprintf(stdout, "Description: one progress per iteration, Stream benchmark\n");
	fprintf(stdout, "Array size:          %zu (elements).\n", array_size);
	fprintf(stdout, "Memory per array:    %.1f MiB.\n",
		(double) memory_size /1024.0/1024.0);
	fprintf(stdout, "Kernel was executed: %ld times.\n", times);
	fprintf(stdout, "Number of threads:   %d\n", num_threads);
	fprintf(stdout, "Progress Time (ns):   %" PRId64 "\n", progress_time);

	for(size_t i = 0; i < 4; i++) {
	fprintf(stdout, "%s Time (s): avg: %11.6f min: %11.6f max: %11.6f\n",
		names[i], 1.0E-09 * sumtime[i]/times, 1.0E-09 * mintime[i],
		1.0E-09 * maxtime[i]);
	fprintf(stdout, "%s Perf (MiB/s): avg: %12.6f best: %12.6f\n", names[i],
		(bytes[i] * 1.0E-06 * memory_size)/ (1.0E-09 * sumtime[i]/times),
		(bytes[i] * 1.0E-06 * memory_size)/ (1.0E-09 * mintime[i]));
	}

#ifdef ENABLE_POST_VALIDATION
	/* validate the benchmark: minimum about of bits should be different. */
	err = 0;
	double ai = 1.0, bi = 2.0, ci = 0.0;
	for(long int i = 0; i < times+1; i++) {
		ci = ai;
		bi = scalar*ci;
		ci = ai+bi;
		ai = bi+scalar*ci;
	}
	for(size_t i = 0; i < array_size && err == 0; i++) {
		err = err || !nrmb_check_double(ai, a[i], 2);
		err = err || !nrmb_check_double(bi, b[i], 2);
		err = err || !nrmb_check_double(ci, c[i], 2);
	}

	if(err)
		fprintf(stdout, "VALIDATION FAILED!!!!\n");
	else
		fprintf(stdout, "VALIDATION PASSED!!!!\n");
	return err;
#else
	fprintf(stdout, "VALIDATION disabled\n");
	return 0;
#endif
}
