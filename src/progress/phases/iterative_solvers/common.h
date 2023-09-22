/*******************************************************************************
 * Copyright 2021 UChicago Argonne, LLC.
 * (c.f. AUTHORS, LICENSE)
 *
 * This file is part of the NRM Benchmarks project.
 * For more info, see https://github.com/anlsys/nrm-benchmarks
 *
 * SPDX-License-Identifier: BSD-3-Clause
 ******************************************************************************/

void initialize_symmetric_positive_good_conditioning(double *A, double *b, double *x, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            double val = 3.2;
            A[i * n + j] = val;
            A[j * n + i] = val;
        }
        A[i * n + i] += n; // Ensure diagonal dominance, hence positive-definiteness
    }

    for (int i = 0; i < n; i++)
    {
        b[i] = 6.5;
    }

    for (int i = 0; i < n; i++)
    {
        x[i] = 0.0;
    }
}

void initialize_symmetric_positive_poor_conditioning(double *A, double *b, double *x, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                A[i * n + j] = (double)(i + 1); // Diagonal elements range from 1 to n
            }
            else
            {
                A[i * n + j] = 0.0; // Off-diagonal elements are zero
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        b[i] = 6.5;
    }

    for (int i = 0; i < n; i++)
    {
        x[i] = 0.0;
    }
}
