/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

// #include "config.h"
// #include "nrm-benchmarks.h"
// #include <nrm.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include <cblas.h>
#include <lapacke.h>

#include "common.h"

static double *A, *b, *x;
int LOG;

// static struct nrm_context *context;

void conjugate_gradient(double *A, double *b, double *x, int n)
{
    double convergence_criteria = 1e-25;
    int total_iterations = 0;

	double *r, *p, *Ap;
    r = (double *)malloc(n * sizeof(double));
    p = (double *)malloc(n * sizeof(double));
    Ap = (double *)malloc(n * sizeof(double));

    cblas_dcopy(n, b, 1, r, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, -1.0, A, n, x, 1, 1.0, r, 1);

    cblas_dcopy(n, r, 1, p, 1);

    double old_residual, residual;
    old_residual = cblas_ddot(n, r, 1, r, 1);

    for (int iter = 0; iter <= n; iter++)
    {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, p, 1, 0.0, Ap, 1);

        double alpha = old_residual / cblas_ddot(n, p, 1, Ap, 1);

#pragma omp parallel for
        for (int j = 0; j < n; j++)
        {
            x[j] = x[j] + alpha * p[j];
            r[j] = r[j] - alpha * Ap[j];
        }

        residual = cblas_ddot(n, r, 1, r, 1);
        if (sqrt(residual) < convergence_criteria)
            break;

#pragma omp parallel for
        for (int j = 0; j < n; j++)
        {
            p[j] = r[j] + (residual / old_residual) * p[j];
        }

        old_residual = residual;
        if (LOG)
        {
            printf("Step: %d Error: %.11lf\n", iter, sqrt(old_residual));
        }
		total_iterations = iter;
    }

    free(r);
    free(p);
    free(Ap);
	
	printf("CG total iterations: %d\n", total_iterations);
}

int main(int argc, char *argv[])
{
    assert(argc == 4);

    int n = atoi(argv[1]);
    char *conditionning = argv[2];
    LOG = atoi(argv[3]);

    A = (double *)malloc(n * n * sizeof(double));
    b = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));

    if (strcmp(conditionning, "good") == 0)
    {
        initialize_symmetric_positive_good_conditioning(A, b, x, n);
    }
    else if (strcmp(conditionning, "poor") == 0)
    {
        initialize_symmetric_positive_poor_conditioning(A, b, x, n);
    }
    else
    {
        fprintf(stderr, "Unknown conditioning option\n");
        exit(EXIT_FAILURE);
    }

	struct timeval start, finish;
	double time;

	gettimeofday(&start, NULL);
    conjugate_gradient(A, b, x, n);
	gettimeofday(&finish, NULL);

    time = (finish.tv_sec - start.tv_sec);
    time += (finish.tv_usec - start.tv_usec)/1e6;
	printf("CG time: %f\n", time);

	free(A);
    free(b);
    free(x);

    return 0;
}
