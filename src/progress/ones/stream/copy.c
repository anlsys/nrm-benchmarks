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

static double *a, *b;
static struct nrm_context *context;

int main(int argc, char **argv)
{
	/* configuration parameters:
	 * - array size in number of double elements
	 * - number of times to run through the benchmark
	 */
	size_t array_size;
	long int times;

	/* needed for performance measurement */
	int64_t sumtime = 0, mintime = INT64_MAX, maxtime = 0;
	nrmb_time_t start, end;
	size_t memory_size;
	int num_threads;

	/* retrieve the size of the allocation and the number of time
	 * to loop through the kernel.
	 */
	assert(argc == 3);
	array_size = strtoull(argv[1], NULL, 0);
	assert(!errno);
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

#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
	{
		a[i] = 1.0;
		b[i] = 2.0;
	}

	/* NRM Context init */
	context = nrm_ctxt_create();
	nrm_init(context, argv[0], 0, 0);

	/* one run of the benchmark for free, warms up the memory */
#pragma omp parallel for
	for(size_t i = 0; i < array_size; i++)
		b[i] = a[i];

	/* this version of the benchmarks reports one progress each time it goes
	 * through the entire array.
	 */
	nrm_send_progress(context, 1);


	for(long int iter = 0; iter < times; iter++)
	{
		int64_t time;
		nrmb_gettime(&start);

		/* the actual benchmark */
#pragma omp parallel for
		for(size_t i = 0; i < array_size; i++)
			b[i] = a[i];

		nrmb_gettime(&end);
		nrm_send_progress(context, 1);

		time = nrmb_timediff(&start, &end);
		sumtime += time;
		mintime = NRMB_MIN(time, mintime);
		maxtime = NRMB_MAX(time, maxtime);
	}

	nrm_fini(context);
	nrm_ctxt_delete(context);

	/* compute stats */

	/* report the configuration and timings */
	fprintf(stdout, "NRM Benchmarks:      %s\n", argv[0]);
	fprintf(stdout, "Version:             %s\n", PACKAGE_VERSION);
	fprintf(stdout, "Description: one progress per iteration, Copy benchmark\n");
	fprintf(stdout, "Array size:          %zu (elements).\n", array_size);
	fprintf(stdout, "Memory per array:    %.1f MiB.\n",
		(double) memory_size /1024.0/1024.0);
	fprintf(stdout, "Kernel was executed: %ld times.\n", times);
	fprintf(stdout, "Number of threads:   %d\n", num_threads);
	fprintf(stdout, "Time (s): avg:       %11.6f min:  %11.6f max: %11.6f\n",
		1.0E-09 * sumtime/times, 1.0E-09 * mintime, 1.0E-09 * maxtime);
	fprintf(stdout, "Perf (MiB/s): avg:   %12.6f best: %12.6f\n",
		(2.0E-06 * memory_size)/ (1.0E-09 * sumtime/times),
		(2.0E-06 * memory_size)/ (1.0E-09 * mintime));

#ifdef ENABLE_POST_VALIDATION
	/* validate the benchmark: for a copy, the minimum about of bits should
	 * be different.
	 */
	err = 0;
	for(size_t i = 0; i < array_size && err == 0; i++)
		err = err || !nrmb_check_double(a[i], b[i], 2);

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
