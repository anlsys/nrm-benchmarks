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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>

#include <cblas.h>

#include "common.h"

static double *A, *b, *x;
int LOG;

static struct nrm_context *context;

void bicgstab(double *A, double *b, double *x, int n, double criteria)
{
    int total_iterations = 0;

    double *r = (double *)malloc(n * sizeof(double));
    double *r_hat = (double *)malloc(n * sizeof(double));
    double *v = (double *)calloc(n, sizeof(double));
    double *p = (double *)malloc(n * sizeof(double));
    double *s = (double *)malloc(n * sizeof(double));
    double *t = (double *)malloc(n * sizeof(double));
    double alpha, omega, rho, rho_prime = 1.0;

    nrm_send_progress(context, 1);

    // Initial residual
    cblas_dcopy(n, b, 1, r, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, -1.0, A, n, x, 1, 1.0, r, 1);
    cblas_dcopy(n, r, 1, r_hat, 1);

    nrm_send_progress(context, 1);

    for (int iter = 0; iter < n; ++iter)
    {
        rho = cblas_ddot(n, r_hat, 1, r, 1);
        if (fabs(rho) < criteria)
            break;

        if (iter == 0)
            cblas_dcopy(n, r, 1, p, 1);
        else
        {
            double beta = (rho / rho_prime) * (alpha / omega);

#pragma omp parallel for
            for (int i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, p, 1, 0.0, v, 1);

        alpha = rho / cblas_ddot(n, r_hat, 1, v, 1);

#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            s[i] = r[i] - alpha * v[i];
        }

        double s_norm = cblas_dnrm2(n, s, 1);
        if (s_norm < 1e-10)
        {
            cblas_daxpy(n, alpha, p, 1, x, 1);
            break;
        }

        cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, s, 1, 0.0, t, 1);

        omega = cblas_ddot(n, t, 1, s, 1) / cblas_ddot(n, t, 1, t, 1);

#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            x[i] = x[i] + alpha * p[i] + omega * s[i];
        }

#pragma omp parallel for
        for (int i = 0; i < n; i++)
        {
            r[i] = s[i] - omega * t[i];
        }

        rho_prime = rho;
        if (LOG)
        {
            printf("Step: %d Error: %.11lf\n", iter, fabs(rho_prime));
        }
		total_iterations = iter;
	nrm_send_progress(context, 1);
    }

    free(r);
    free(r_hat);
    free(v);
    free(p);
    free(s);
    free(t);

	printf("BiCGStab total iterations: %d\n", total_iterations);
}

int main(int argc, char *argv[])
{
    assert(argc == 5);

    int n = atoi(argv[1]);
    char *conditionning = argv[2];
    LOG = atoi(argv[3]);
    int criteria = atoi(argv[4]);

    A = (double *)malloc(n * n * sizeof(double));
    b = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));

    double convergence = pow(1, -criteria);

    context = nrm_ctxt_create();
    nrm_init(context, argv[0], 0, 0);
    nrm_send_progress(context, 1);

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
    bicgstab(A, b, x, n, convergence);
    gettimeofday(&finish, NULL);

    nrm_fini(context);
    nrm_ctxt_delete(context);

    time = (finish.tv_sec - start.tv_sec);
    time += (finish.tv_usec - start.tv_usec)/1e6;
	printf("BiCGStab time: %f\n", time);
    
	free(A);
    free(b);
    free(x);

    return 0;
}
