
/*******************************************************************************
 *
 * File update.c
 *
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to perform gauge field sweeps.
 *
 * Externally accessible functions:
 *
 * void update(int iup,int nup,int utype,int stype)
 *      Performs nup sweeps of local updates of the type defined by utype.
 *      stype defines the sweep type.
 *
 * utype==1 : 1 Metropolis update
 *
 * stype==0 : Sweeps with all links in order of appearance
 * stype==1 : Sweeps with randomly chosen points
 *
 *******************************************************************************/

#define UPDATE_C

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ranlxd.h"
#include "headers.h"
#include "modules.h"

#if (MASTER_FIELD == 1)
void gaugefieldUpdate(int iup, int nup, int utype, int stype)
{
    int n, in, dir, idir, nn;
    double nsum, msum, r[2];
    double iaccRate, t1, t2;

    checkpoint("update -- in");

    error(utype != 0 && utype != 1, "update [update.c]", "Unknown update type! Only "
                                                         "Metropolis (SIM_TYPE=0) and heatbath (SIM_TYPE=1) implemented!");
    error(stype != 0 && stype != 1, "update [update.c]", "Unknown sweep type! "
                                                         "Only sequential (SWEEP_TYPE=0) and random (SWEEP_TYPE=1) "
                                                         "sweeps implemented!");

    t1 = getTime();

    msum = 0.;
    for (nn = 0; nn < nup; nn++)
    {
        nsum = 0;
        for (n = 0; n < VOL; n++)
        {
            if (stype == 1)
            {
                ranlxd(r, 1);
                in = (int)(VOL * r[0]);
                if (in == VOL)
                    in = VOL - 1;
            }
            else
                in = n;
            
            in = i[in]; //only difference to non master-field sim
            for (dir = 0; dir < DIM; dir++)
            {
                if (stype == 1)
                {
                    idir = (int)(DIM * r[1]);
                    if (idir == DIM)
                        idir = DIM - 1;
                }
                else
                    idir = dir;

                if (utype == 0)
                    nsum += localMetropolisUpdate(in, dir, 1);
                else if (utype == 1)
                    nsum += localHeatbathUpdate(in, dir, 1);
            }
        }

        msum += nsum / (VOL * DIM);
    }
    iaccRate = msum / nup;
    project_gfield_to_sun(pu); //?????

    t2 = getTime();

    logging("\nUpdate %d done:\n"
            "(algorithm %d; sweep %d; acc=%.3f)\n"
            "update took %.9f seconds\n",
            iup, utype, stype, iaccRate, t2 - t1);

    checkpoint("update -- out");
}
#elif (MASTER_FIELD == 0)
void gaugefieldUpdate(int iup, int nup, int utype, int stype)
{
    int n, in, dir, idir, nn;
    double nsum, msum, r[2];
    double iaccRate, t1, t2;

    checkpoint("update -- in");

    error(utype != 0 && utype != 1, "update [update.c]", "Unknown update type! Only "
                                                         "Metropolis (SIM_TYPE=0) and heatbath (SIM_TYPE=1) implemented!");
    error(stype != 0 && stype != 1, "update [update.c]", "Unknown sweep type! "
                                                         "Only sequential (SWEEP_TYPE=0) and random (SWEEP_TYPE=1) "
                                                         "sweeps implemented!");

    t1 = getTime();

    msum = 0.;
    for (nn = 0; nn < nup; nn++)
    {
        nsum = 0;
        for (n = 0; n < VOL; n++)
        {
            if (stype == 1)
            {
                ranlxd(r, 1);
                in = (int)(VOL * r[0]);
                if (in == VOL)
                    in = VOL - 1;
            }
            else
                in = n;
            for (dir = 0; dir < DIM; dir++)
            {
                if (stype == 1)
                {
                    idir = (int)(DIM * r[1]);
                    if (idir == DIM)
                        idir = DIM - 1;
                }
                else
                    idir = dir;

                if (utype == 0)
                    nsum += localMetropolisUpdate(in, dir, 1);
                else if (utype == 1)
                    nsum += localHeatbathUpdate(in, dir, 1);
            }
        }

        msum += nsum / (VOL * DIM);
    }
    iaccRate = msum / nup;
    project_gfield_to_sun(pu);

    t2 = getTime();

    logging("\nUpdate %d done:\n"
            "(algorithm %d; sweep %d; acc=%.3f)\n"
            "update took %.9f seconds\n",
            iup, utype, stype, iaccRate, t2 - t1);

    checkpoint("update -- out");
}
#endif