
/*******************************************************************************
 *
 * File heatbath.c
 *
 * Copyright (C) 2022 Marco Stilger
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to perform heatbath updates for pure gauge theory in SU(2).
 *
 * Externally accessible functions:
 *
 * void staples(int n,int dir,sun_mat *stap)
 *      Computes the staples for the link starting at point n in direction
 *      dir. The staples matrix is returned via the pointer stap.
 *
 * double localHeatbathUpdate(int n,int dir,int m)
 *        Performs m local Heatbath update steps for the link U_dir(n).
 *        On return it hands back the fraction of successful updates.
 *        For now only works for SUN=2.
 *
 *******************************************************************************/

#define HEATHBATH_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ranlxd.h"
#include "headers.h"
#include "modules.h"

/*
void staples(int n, int dir, sun_mat *stap)
{
    int kk, n1, n2;
    sun_mat un[4];

    sun_zero(un[3]);

    for (kk = 0; kk < DIM; kk++)
    {
        if (kk != dir)
        {
            n1 = neib[n][dir];
            n2 = neib[n][kk];

            sun_mul_dag(un[1], *pu[n1][kk], *pu[n2][dir]);
            sun_mul_dag(un[0], un[1], *pu[n][kk]);
            sun_add(un[2], un[3], un[0]);

            n2 = neib[n][kk + DIM];
            n1 = neib[n2][dir];

            sun_dag(un[0], *pu[n1][kk]);
            sun_mul_dag(un[1], un[0], *pu[n2][dir]);
            sun_mul(un[0], un[1], *pu[n2][kk]);
            sun_add(un[3], un[2], un[0]);
        }
    }

    *stap = un[3];
}
*/

double localHeatbathUpdate(int n, int dir, int m)
{
#if (SUN == 2)
    int im, count = 0;
    double sqrt_det, ud[4], lambda_sqr, angles[2], length, x[4];

    su2mat A, X;
    staples(n, dir, &A);
    su2_det_sqrt(sqrt_det, A);
    if (sqrt_det <= __DBL_EPSILON__)
    {
        /*su2RandomMatrix(pu[n][dir]);*/
        return 0.0;
    }
    su2_dble_div(A, sqrt_det);

    for (im = 0; im < m; im++)
    {
        do
        {
            ranlxd(ud, 4);
            lambda_sqr = -1.0 / (2.0 * sqrt_det * runParams.beta) * (log(1. - ud[1]) + pow(cos(2 * M_PI * (1. - ud[2])), 2) * log(1. - ud[3]));
            count++;
        } while (pow(ud[3], 2) > 1 - pow(lambda_sqr, 2));

        ranlxd(angles, 2);
        angles[0] *= 2.;
        angles[0] -= 1.;
        angles[0] = acos(angles[1]);
        angles[1] *= 2. * PI;

        x[0] = 1 - 2.0*lambda_sqr;
        length = sqrt(1 - pow(x[0], 2));
        x[1] = length * sin(angles[0]) * cos(angles[1]);
        x[2] = length * sin(angles[0]) * sin(angles[1]);
        x[3] = length * cos(angles[0]);

        mk_dble_array_su2(x, X);

#if DEBUG == 1
        double det_1, det_2;
        su2mat X_dag, res;
        su2_dag(X_dag, X);
        su2_mat_mul(res, X_dag, X);
        su2_det(det_2, res);
        su2_det(det_1, X);
        printf("%3.20f   %3.20f\n", det_1, det_2);
#endif
        su2_mat_mul_dag(*pu[n][dir], X, A);
    }
    return (double)m / count;

#elif (SUN == 3)

#endif
}