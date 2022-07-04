
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

double generateSU2InHeatbath(su2mat *X, double alpha)
{
    int count = 0;
    double ud[4], lambda_sqr, angles[2], length, x[4];

    do
    {
        ranlxd(ud, 4);
        lambda_sqr = -1.0 / alpha * (log(1. - ud[1]) + pow(cos(2 * M_PI * (1. - ud[2])), 2) * log(1. - ud[3]));
        count++;
    } while (pow(ud[3], 2) > 1 - 0.5 * pow(lambda_sqr, 2));

    ranlxd(angles, 2);
    angles[0] *= 2.;
    angles[0] -= 1.;
    angles[0] = acos(angles[1]);
    angles[1] *= 2. * PI;

    x[0] = 1 - lambda_sqr;
    length = sqrt(1 - pow(x[0], 2));
    x[1] = length * sin(angles[0]) * cos(angles[1]);
    x[2] = length * sin(angles[0]) * sin(angles[1]);
    x[3] = length * cos(angles[0]);

    mk_dble_array_su2(x, *X);

#if DEBUG == 1
    double det_1, det_2;
    su2mat X_dag, res;
    su2_dag(X_dag, X);
    su2_mat_mul(res, X_dag, X);
    su2_det(det_2, res);
    su2_det(det_1, X);
    printf("%3.20f   %3.20f\n", det_1, det_2);
#endif
    return count;
}

double localHeatbathUpdate(int n, int dir, int m)
{

    int count = 0;
    double sqrt_det, alpha;

#if (SUN == 2)
    su2mat A, X, U;

    staples(n, dir, &A);
    su2_det_sqrt(sqrt_det, A);

    if (sqrt_det <= __DBL_EPSILON__)
    {
        /*su2RandomMatrix(pu[n][dir]);*/
        return 0.0;
    }

    su2_dble_div(A, sqrt_det);

    alpha = sqrt_det * runParams.beta;

    count += generateSU2InHeatbath(&X, alpha);
    su2_mat_mul_dag(U, X, A); /*U = X*V^dag*/
    *pu[n][dir] = U;
    return 1.0 / count;

#elif (SUN == 3)

    su3mat heatbath_su3, old_link, staple_sum, Temp_1, Temp_2;
    su2mat heatbath_su2, generated, new_link_su2;

    staples(n, dir, &staple_sum);
    old_link = *pu[n][dir];

    for (int i = 1; i <= 3; i++)
    {
        su3_mat_mul(heatbath_su3, old_link, staple_sum);
        switch (i)
        {
        case 1:
            extr_sub1(heatbath_su2, heatbath_su3);
            break;
        case 2:
            extr_sub2(heatbath_su2, heatbath_su3);
            break;
        case 3:
            extr_sub3(heatbath_su2, heatbath_su3);
            break;
        }
        su2_det_sqrt(sqrt_det, heatbath_su2);

        if (sqrt_det <= __DBL_EPSILON__)
        {
            /*su3RandomMatrix(pu[n][dir]);*/
            return 0.0;
        }
        su2_dble_div(heatbath_su2, sqrt_det);

        alpha = sqrt_det * runParams.beta * 2.0 / 3.0;

        count += generateSU2InHeatbath(&generated, alpha);
        su2_mat_mul_dag(new_link_su2, generated, heatbath_su2);

        switch (i)
        {
        case 1:
            create_su3_sub1(Temp_1, new_link_su2);
            break;
        case 2:
            create_su3_sub2(Temp_1, new_link_su2);
            break;
        case 3:
            create_su3_sub3(Temp_1, new_link_su2);
            break;
        }

        su3_mat_mul(Temp_2, Temp_1, heatbath_su3);
        heatbath_su3 = Temp_2;
        su3_mat_mul(Temp_2, Temp_1, old_link);
        old_link = Temp_2;
    }

    *pu[n][dir] = old_link;

    return 3.0 / count;
#endif
}