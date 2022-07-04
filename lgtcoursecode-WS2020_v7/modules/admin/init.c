/*******************************************************************************
 *
 * File init.c
 *
 * Copyright (C) 2020 Alessandro Sciarra
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines to initialise the fields and the error checks.
 *
 * Externally accessible functions:
 *
 * void initProgram(int itype)
 *      Checks runParameters.
 *      Initialises some constants etc..
 * 
 * void initArrayOfNeighbours(void)
 *      Initialises the neib arrays.
 *
 * void initGaugeField(int flag)
 *      Initialises the gauge field, i.e. the global sun_mat *pu[VOL][DIM]
 *      with either unity or random SU(N) matrices depending on the flag value.
 *
 * void releaseGaugeField(void)
 *      Frees the memory allocated for the gauge field.
 * 
 * void initArrayOfSublattices(void)
 *      For master-field simulation. Initalizes array sl.
 * 
 * void initArrayOfPermutation(int size, int perm[][size])
 *      For master-field simulation. Saves all permutations of integers from 0 to size-1 in perm.
 * 
 * void updateArrayOfIndexMF(int)
 *      For master-field simulation. fills array index_mf[VOL2] depending on sublattice index. 
 *      Has to be called when going to different sublattice.
 * 
 * void initArrayOfI(void)
 *      For master-field simulation. initializes array i.
 * 
 * void initGlobalArrays(void)
 *      For master-field simulation. Wrapper for all init functions. not including updateArrayOfIndexMF.
 *******************************************************************************/

#define INIT_C

#include "modules.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>

#define DEBUG_INIT 0

void initProgram(int itype)
{

    error((runParams.numConfs <= 0) || (runParams.decorSteps < 0) || (runParams.numThermConfs < 0) ||
              ((((runParams.numConfs - runParams.numThermConfs) % runParams.writeConfsFreq) != 0) && (runParams.writeConfsFreq > 0)),
          "initProgram [init.c]", "Error in run parameters");

    PI = 2.0 * asin(1.0);

    initTwoPi();
}

#if (MASTER_FIELD == 1)
void initArrayOfSubLattices(void)
{
    int n, dir, ns = 0;
    int volumeFactor[DIM], latticeExtent[DIM], latticeExtent_mf[DIM], coord[DIM];

    checkpoint("initArrayOfSubLattices -- in");
#if (DIM == 2)
    volumeFactor[0] = LENGS1_MF;
    volumeFactor[1] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    coord[0] = 0;
    coord[1] = 0;
#elif (DIM == 3)
    volumeFactor[0] = LENGS2_MF * LENGS1_MF;
    volumeFactor[1] = LENGS2_MF;
    volumeFactor[2] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
#elif (DIM == 4)
    volumeFactor[0] = LENGS3_MF * LENGS2_MF * LENGS1_MF;
    volumeFactor[1] = LENGS3_MF * LENGS2_MF;
    volumeFactor[2] = LENGS3_MF;
    volumeFactor[3] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    latticeExtent[3] = LENGS3;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
    latticeExtent_mf[3] = LENGS3_MF;
    coord[0] = 0;
    coord[1] = 0;
    coord[2] = 0;
    coord[3] = 0;
#endif

    for (n = 0; n < VOL_SL; n++)
    {
        sl[n] = ns;
        for (dir = DIM - 1; dir >= 0; dir--)
        {
            if (coord[dir] + latticeExtent[dir] < latticeExtent_mf[dir])
            {
                coord[dir] += latticeExtent[dir];
                break;
            }
            else
            {
                coord[dir] = 0;
            }
        }
        ns = 0;
        for (dir = 0; dir < DIM; dir++)
        {
            ns += volumeFactor[dir] * coord[dir];
        }
    }

#if (DEBUG_INIT == 1)
    logging("\nix");
    for (n = 0; n < VOL_SL; n++)
    {
        logging("\n%d", n);
        logging("\t%d", sl[n]);
    }
    logging("\n");
#endif

    checkpoint("initArrayOfSubLattices -- out");
    return;
}

/************************************************************************************************************************/

static void initArrayOfPermutations(int size, int perm[][size])
{
    int init[size], p[size + 1];
    int num_perms = fact(size);
    int i, j, k, count;
    for (i = 0; i < size; i++)
    {
        init[i] = i;
    }
    for (i = 0; i < size; i++)
    {
        p[i] = 0;
    }

    count = 0;
    for (i = 0; i < size; i++)
    {
        perm[count][i] = init[i];
    }

    i = 1;
    while (i < size)
    {
        if (p[i] < i)
        {
            if (i % 2 == 0)
            {
                swap(&init[0], &init[i]);
            }
            else
            {
                swap(&init[p[i]], &init[i]);
            }

            count++;
            for (j = 0; j < num_perms; j++)
            {
                for (k = 0; k < size; k++)
                {
                    perm[count][k] = init[k];
                }
            }

            p[i]++;
            i = 1;
        }
        else
        {
            p[i] = 0;
            i++;
        }
    }
}

void updateArrayOfIndexMF(int n_sl)
{
    long int latticeExtent[DIM], latticeExtent_mf[DIM];
    long int volumeFactor[DIM], volumeFactor_big[DIM], volumeFactor_mf[DIM];
    long int ix[DIM], ix_big[DIM], ix_mf[DIM], coord_big[DIM], coord_mf[DIM];
    long int in, dir, dir2, index, index_m, border_count, num_perms;
    bool border[2 * DIM];
    double t1, t2;

    checkpoint("initArrayOfIndexMF -- in");
    t1 = getTime();

#if (DIM == 2)
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    volumeFactor_mf[0] = LENGS1_MF;
    volumeFactor_mf[1] = 1;
    volumeFactor_big[0] = (LENGS1 + 2);
    volumeFactor_big[1] = 1;
    volumeFactor[0] = LENGS1;
    volumeFactor[1] = 1;
#elif (DIM == 3)
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
    volumeFactor_mf[0] = LENGS1_MF * LENGS2_MF;
    volumeFactor_mf[1] = LENGS2_MF;
    volumeFactor_mf[2] = 1;
    volumeFactor_big[0] = (LENGS1 + 2) * (LENGS2 + 2);
    volumeFactor_big[1] = (LENGS2 + 2);
    volumeFactor_big[2] = 1;
    volumeFactor[0] = LENGS1 * LENGS2;
    volumeFactor[1] = LENGS2;
    volumeFactor[2] = 1;
#elif (DIM == 4)
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    latticeExtent[3] = LENGS3;
    latticeExtent_mf[0] = LENGT_MF;
    latticeExtent_mf[1] = LENGS1_MF;
    latticeExtent_mf[2] = LENGS2_MF;
    latticeExtent_mf[3] = LENGS3_MF;
    volumeFactor_mf[0] = LENGS1_MF * LENGS2_MF * LENGS3_MF;
    volumeFactor_mf[1] = LENGS2_MF * LENGS3_MF;
    volumeFactor_mf[2] = LENGS3_MF;
    volumeFactor_mf[3] = 1;
    volumeFactor_big[0] = (LENGS1 + 2) * (LENGS2 + 2) * (LENGS3 + 2);
    volumeFactor_big[1] = (LENGS2 + 2) * (LENGS3 + 2);
    volumeFactor_big[2] = (LENGS3 + 2);
    volumeFactor_big[3] = 1;
    volumeFactor[0] = LENGS1 * LENGS2 * LENGS3;
    volumeFactor[1] = LENGS2 * LENGS3;
    volumeFactor[2] = LENGS3;
    volumeFactor[3] = 1;
#endif
    for (dir = 0; dir < DIM; dir++)
    {
        ix[dir] = 0;
        ix_big[dir] = 0;
        ix_mf[dir] = 0;
        coord_big[dir] = 0;
        coord_mf[dir] = 0;
    }

    index = 0;
    for (in = 0; in < VOL; in++) /* calculate index_mf for core part */
    {
        index_mf[i[in]] = index + sl[n_sl];

        for (dir = DIM - 1; dir >= 0; dir--) /* calculate coordinates of small lattice*/
        {
            if (ix[dir] + 1 != latticeExtent[dir])
            {
                ix[dir]++;
                break;
            }
            else
            {
                ix[dir] = 0;
            }
        }
        index = 0;
        for (int dir = 0; dir < DIM; dir++) /* calculate lattice site index from coordinate */
        {
            index += ix[dir] * volumeFactor_mf[dir];
        }
    }

    for (in = 0; in < VOL; in++) /* calculate index_mf for outer layers*/
    {
        index = in;
        for (dir = 0; dir < DIM; dir++)
        {
            ix[dir] = index / volumeFactor[dir];
            index = index % volumeFactor[dir];
        }

        index = i[in];
        for (dir = 0; dir < DIM; dir++)
        {
            ix_big[dir] = index / volumeFactor_big[dir];
            index = index % volumeFactor_big[dir];
        }

        index = index_mf[i[in]];
        for (dir = 0; dir < DIM; dir++)
        {
            ix_mf[dir] = index / volumeFactor_mf[dir];
            index = index % volumeFactor_mf[dir];
        }

        border_count = 0;
        for (dir = 0; dir < DIM; dir++)
        {
            border[dir] = false;
            border[dir + DIM] = false;
            if (ix[dir] == 0)
            {
                border[dir + DIM] = true;
                border_count++;
            }
            else if (ix[dir] == latticeExtent[dir] - 1)
            {
                border[dir] = true;
                border_count++;
            }
        }

        num_perms = fact(border_count);
        int borders[border_count];
        int perms[num_perms][border_count];
        border_count = 0;
        for (dir = 0; dir < DIM; dir++)
        {
            if (border[dir + DIM] == true)
            {
                borders[border_count] = (dir + DIM);
                border_count++;
            }
            else if (border[dir] == true)
            {
                borders[border_count] = dir;
                border_count++;
            }
        }

        if (border_count == 1)
        {
            for (int i = 0; i < num_perms; i++)
            {
                for (int j = 0; j < border_count; j++)
                {
                    perms[i][j] = perm1[i][j];
                }
            }
        }
        else if (border_count == 2)
        {
            for (int i = 0; i < num_perms; i++)
            {
                for (int j = 0; j < border_count; j++)
                {
                    perms[i][j] = perm2[i][j];
                }
            }
        }
        else if (border_count == 3)
        {
            for (int i = 0; i < num_perms; i++)
            {
                for (int j = 0; j < border_count; j++)
                {
                    perms[i][j] = perm3[i][j];
                }
            }
        }
        else if (border_count == 4)
        {
            for (int i = 0; i < num_perms; i++)
            {
                for (int j = 0; j < border_count; j++)
                {
                    perms[i][j] = perm4[i][j];
                }
            }
        }
        //initArrayOfPermutations(border_count, perms);

        for (dir = 0; dir < DIM; dir++)
        {
            coord_mf[dir] = ix_mf[dir];
            coord_big[dir] = ix_big[dir];
        }

        for (int i = 0; i < num_perms; i++)
        {
            for (int j = 0; j < border_count; j++)
            {
                dir = borders[perms[i][j]];
                if (dir < DIM)
                {
                    ix_mf[dir]++;
                    ix_big[dir]++;
                }
                else
                {
                    dir -= DIM;
                    ix_mf[dir]--;
                    ix_big[dir]--;
                }
                ix_mf[dir] = (ix_mf[dir] + latticeExtent_mf[dir]) % latticeExtent_mf[dir];

                index_m = 0;
                index = 0;
                for (dir2 = 0; dir2 < DIM; dir2++)
                {
                    index += volumeFactor_big[dir2] * ix_big[dir2];
                    index_m += volumeFactor_mf[dir2] * ix_mf[dir2];
                }
                index_mf[index] = index_m;
            }
            for (dir = 0; dir < DIM; dir++)
            {
                ix_mf[dir] = coord_mf[dir];
                ix_big[dir] = coord_big[dir];
            }
        }
    }

#if (DEBUG_INIT == 1)
    logging("\n\n\nin\t-->\ti\n\n");

    for (in = 0; in < VOL2; in++)
    {
        logging("%d\t-->\t%d\n", in, index_mf[in]);
    }

#endif

    t2 = getTime();
    logging("\nArray index_mf updated at position %d (took %.3e sec)!\n", n_sl, t2 - t1);
    checkpoint("initArrayOfIndexMF -- out");
    return;
}

void initArrayOfNeighbours(void) // works with VOL2 stuff
{
    int i, index, dir1, dir2;
    int coord[DIM], newCoord[DIM];
    int volumeFactor[DIM];
    int latticeExtent[DIM], volumeOtherDirs[DIM];

    checkpoint("initArrayOfNeighbours -- in");

#if (DIM == 2)
    volumeFactor[0] = (LENGS1 + 2);
    volumeFactor[1] = 1;
    latticeExtent[0] = (LENGT + 2);
    latticeExtent[1] = (LENGS1 + 2);
    volumeOtherDirs[0] = SVOL2;
    volumeOtherDirs[1] = (LENGT + 2);
#elif (DIM == 3)
    volumeFactor[0] = (LENGS1 + 2) * (LENGS2 + 2);
    volumeFactor[1] = (LENGS2 + 2);
    volumeFactor[2] = 1;
    latticeExtent[0] = (LENGT + 2);
    latticeExtent[1] = (LENGS1 + 2);
    latticeExtent[2] = (LENGS2 + 2);
    volumeOtherDirs[0] = SVOL2;
    volumeOtherDirs[1] = (LENGT + 2) * (LENGS2 + 2);
    volumeOtherDirs[2] = (LENGT + 2) * (LENGS1 + 2);
#elif (DIM == 4)
    volumeFactor[0] = (LENGS1 + 2) * (LENGS2 + 2) * (LENGS3 + 2);
    volumeFactor[1] = (LENGS2 + 2) * (LENGS3 + 2);
    volumeFactor[2] = (LENGS3 + 2);
    volumeFactor[3] = 1;
    latticeExtent[0] = (LENGT + 2);
    latticeExtent[1] = (LENGS1 + 2);
    latticeExtent[2] = (LENGS2 + 2);
    latticeExtent[3] = (LENGS3 + 2);
    volumeOtherDirs[0] = SVOL2;
    volumeOtherDirs[1] = (LENGT + 2) * (LENGS2 + 2) * (LENGS3 + 2);
    volumeOtherDirs[2] = (LENGT + 2) * (LENGS1 + 2) * (LENGS3 + 2);
    volumeOtherDirs[3] = (LENGT + 2) * (LENGS1 + 2) * (LENGS2 + 2);
#endif

    for (dir1 = 0; dir1 < DIM; dir1++)
    {
        latParams.linearExtent[dir1] = latticeExtent[dir1];
        latParams.volumeOtherDirs[dir1] = volumeOtherDirs[dir1];
        latParams.volumeFactor[dir1] = volumeFactor[dir1];
    }

    for (i = 0; i < VOL2; i++)
    {
        // Compute coordinates in all directions for site index i
        for (dir1 = 0, index = i; dir1 < DIM; dir1++)
        {
            coord[dir1] = index / volumeFactor[dir1];
            index = index % volumeFactor[dir1];
        }
        for (dir1 = 0; dir1 < DIM; dir1++)
        {
            for (dir2 = 0; dir2 < DIM; dir2++)
                newCoord[dir2] = coord[dir2];

            // Compute coordinates of nearest neighbour in positive directions
            newCoord[dir1] = (coord[dir1] + 1) % latticeExtent[dir1];
            for (dir2 = 0, neib[i][dir1] = 0; dir2 < DIM; dir2++)
                neib[i][dir1] += newCoord[dir2] * volumeFactor[dir2];

            // Compute coordinates of nearest neighbour in negative directions
            newCoord[dir1] = (coord[dir1] + latticeExtent[dir1] - 1) % latticeExtent[dir1];
            for (dir2 = 0, neib[i][dir1 + DIM] = 0; dir2 < DIM; dir2++)
                neib[i][dir1 + DIM] += newCoord[dir2] * volumeFactor[dir2];
        }
    }

#if (DEBUG_INIT == 1)
    logging("\nix");
    for (dir1 = 0; dir1 < DIM; dir1++)
        logging("\tn%d\tn%d", dir1, dir1 + DIM);
    for (i = 0; i < VOL2; i++)
    {
        logging("\n%d", i);
        for (dir1 = 0; dir1 < DIM; dir1++)
            logging("\t%d\t%d", neib[i][dir1], neib[i][dir1 + DIM]);
    }
    logging("\n");
#endif

    checkpoint("initArrayOfNeighbours -- out");
    return;
}

void initArrayOfI(void) // works with VOL2 stuff
{
    checkpoint("initArrayOfI -- in");

    int n, count, dir;
    bool border;

    count = 0;
    for (n = 0; n < VOL2; n++)
    {
        border = false;
        for (dir = 0; dir < DIM; dir++)
        {
            if (neib[n][dir] < n || neib[n][dir + DIM] > n)
                border = true;
        }
        if (!border)
        {
            i[count] = n;
            count++;
        }
    }
#if (DEBUG_INIT == 1)
    logging("n\ti(n)\n");
    for (n = 0; n < VOL; n++)
    {
        logging("%d\t%d\n", n, i[n]);
    }
#endif

    checkpoint("initArrayofI -- out");
    return;
}

void initGlobalArrays(void)
{
    initArrayOfNeighbours();
    initArrayOfSubLattices();
    initArrayOfI();
    initArrayOfPermutations(1, perm1);
    initArrayOfPermutations(2, perm2);
    initArrayOfPermutations(3, perm3);
    initArrayOfPermutations(4, perm4);
    return;
}

void initGaugeField(int flag)
{
    sun_mat *iu = NULL, *zu, ini;
    int ind, dir;

    checkpoint("init_gauge");

    error((flag != 0) && (flag != 1), "init_gauge [start.c]",
          "Wrong flag! Cannot initiate gauge field!");

    iu = malloc(DIM * VOL2 * sizeof(sun_mat));
    error(iu == NULL, "init_gauge [start.c]",
          "Unable to allocate gauge field array!");

    if (flag == 0)
    {
        sun_unit(ini);
    }
    for (ind = 0, zu = iu; ind < VOL2; ind++)
    {
        for (dir = 0; dir < DIM; dir++, zu++)
        {
            pu[ind][dir] = zu;
            if (flag == 0)
            {
                *pu[ind][dir] = ini;
            }
            else if (flag == 1)
            {
#if (SUN == 3)
                su3RandomMatrix(pu[ind][dir]);
#elif (SUN == 2)
                su2RandomMatrix(pu[ind][dir]);
#endif
            }
        }
    }
}

void allocateGaugeField(sun_mat *u[VOL2][DIM])
{
    sun_mat *iu = NULL, *zu;
    int ind, dir;

    checkpoint("allocate_gauge");

    iu = malloc(DIM * VOL2 * sizeof(sun_mat));
    error(iu == NULL, "allocate_gauge [start.c]",
          "Unable to allocate gauge field array!");

    for (ind = 0, zu = iu; ind < VOL2; ind++)
    {
        for (dir = 0; dir < DIM; dir++, zu++)
        {
            u[ind][dir] = zu;
        }
    }
}

void deallocateGaugeField(sun_mat *u[VOL2][DIM])
{
    free(u[0][0]);
}

void copyGaugeField(sun_mat *u1[VOL2][DIM], sun_mat *u2[VOL2][DIM])
{
    int ii, jj;
    for (ii = 0; ii < VOL2; ii++)
        for (jj = 0; jj < DIM; jj++)
            *u2[ii][jj] = *u1[ii][jj];
}

#elif (MASTER_FIELD == 0)
void initArrayOfNeighbours(void)
{
    int i, index, dir1, dir2;
    int coord[DIM], newCoord[DIM];
    int volumeFactor[DIM];
    int latticeExtent[DIM], volumeOtherDirs[DIM];

    checkpoint("initArrayOfNeighbours -- in");

#if (DIM == 2)
    volumeFactor[0] = LENGS1;
    volumeFactor[1] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    volumeOtherDirs[0] = SVOL;
    volumeOtherDirs[1] = LENGT;
#elif (DIM == 3)
    volumeFactor[0] = LENGS1 * LENGS2;
    volumeFactor[1] = LENGS2;
    volumeFactor[2] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    volumeOtherDirs[0] = SVOL;
    volumeOtherDirs[1] = LENGT * LENGS2;
    volumeOtherDirs[2] = LENGT * LENGS1;
#elif (DIM == 4)
    volumeFactor[0] = LENGS1 * LENGS2 * LENGS3;
    volumeFactor[1] = LENGS2 * LENGS3;
    volumeFactor[2] = LENGS3;
    volumeFactor[3] = 1;
    latticeExtent[0] = LENGT;
    latticeExtent[1] = LENGS1;
    latticeExtent[2] = LENGS2;
    latticeExtent[3] = LENGS3;
    volumeOtherDirs[0] = SVOL;
    volumeOtherDirs[1] = LENGT * LENGS2 * LENGS3;
    volumeOtherDirs[2] = LENGT * LENGS1 * LENGS3;
    volumeOtherDirs[3] = LENGT * LENGS1 * LENGS2;
#endif

    for (dir1 = 0; dir1 < DIM; dir1++)
    {
        latParams.linearExtent[dir1] = latticeExtent[dir1];
        latParams.volumeOtherDirs[dir1] = volumeOtherDirs[dir1];
        latParams.volumeFactor[dir1] = volumeFactor[dir1];
    }

    for (i = 0; i < VOL; i++)
    {
        // Compute coordinates in all directions for site index i
        for (dir1 = 0, index = i; dir1 < DIM; dir1++)
        {
            coord[dir1] = index / volumeFactor[dir1];
            index = index % volumeFactor[dir1];
        }
        for (dir1 = 0; dir1 < DIM; dir1++)
        {
            for (dir2 = 0; dir2 < DIM; dir2++)
                newCoord[dir2] = coord[dir2];

            // Compute coordinates of nearest neighbour in positive directions
            newCoord[dir1] = (coord[dir1] + 1) % latticeExtent[dir1];
            for (dir2 = 0, neib[i][dir1] = 0; dir2 < DIM; dir2++)
                neib[i][dir1] += newCoord[dir2] * volumeFactor[dir2];

            // Compute coordinates of nearest neighbour in negative directions
            newCoord[dir1] = (coord[dir1] + latticeExtent[dir1] - 1) % latticeExtent[dir1];
            for (dir2 = 0, neib[i][dir1 + DIM] = 0; dir2 < DIM; dir2++)
                neib[i][dir1 + DIM] += newCoord[dir2] * volumeFactor[dir2];
        }
    }

#if (DEBUG_INIT == 1)
    logging("\nix");
    for (dir1 = 0; dir1 < DIM; dir1++)
        logging("\tn%d\tn%d", dir1, dir1 + DIM);
    for (i = 0; i < VOL; i++)
    {
        logging("\n%d", i);
        for (dir1 = 0; dir1 < DIM; dir1++)
            logging("\t%d\t%d", neib[i][dir1], neib[i][dir1 + DIM]);
    }
    logging("\n");
#endif

    checkpoint("initArrayOfNeighbours -- out");
    return;
}

void initGaugeField(int flag)
{
    sun_mat *iu = NULL, *zu, ini;
    int ind, dir;

    checkpoint("init_gauge");

    error((flag != 0) && (flag != 1), "init_gauge [start.c]",
          "Wrong flag! Cannot initiate gauge field!");

    iu = malloc(DIM * VOL * sizeof(sun_mat));
    error(iu == NULL, "init_gauge [start.c]",
          "Unable to allocate gauge field array!");

    if (flag == 0)
    {
        sun_unit(ini);
    }
    for (ind = 0, zu = iu; ind < VOL; ind++)
    {
        for (dir = 0; dir < DIM; dir++, zu++)
        {
            pu[ind][dir] = zu;
            if (flag == 0)
            {
                *pu[ind][dir] = ini;
            }
            else if (flag == 1)
            {
#if (SUN == 3)
                su3RandomMatrix(pu[ind][dir]);
#elif (SUN == 2)
                su2RandomMatrix(pu[ind][dir]);
#endif
            }
        }
    }
}

void allocateGaugeField(sun_mat *u[VOL][DIM])
{
    sun_mat *iu = NULL, *zu;
    int ind, dir;

    checkpoint("allocate_gauge");

    iu = malloc(DIM * VOL * sizeof(sun_mat));
    error(iu == NULL, "allocate_gauge [start.c]",
          "Unable to allocate gauge field array!");

    for (ind = 0, zu = iu; ind < VOL; ind++)
    {
        for (dir = 0; dir < DIM; dir++, zu++)
        {
            u[ind][dir] = zu;
        }
    }
}

void deallocateGaugeField(sun_mat *u[VOL][DIM])
{
    free(u[0][0]);
}

void copyGaugeField(sun_mat *u1[VOL][DIM], sun_mat *u2[VOL][DIM])
{
    int ii, jj;
    for (ii = 0; ii < VOL; ii++)
        for (jj = 0; jj < DIM; jj++)
            *u2[ii][jj] = *u1[ii][jj];
}

#endif

void allocateFermionField(sun_wferm **f)
{
    sun_wferm *s;
    s = malloc(VOL * sizeof(sun_wferm));
    error(s == NULL, "allocateFermionField [init.c]",
          "Unable to allocate fermion field array!");
    *f = s;
}

void deallocateFermionField(sun_wferm **f)
{
    free(*f);
}

void releaseGaugeField(void)
{
    free(pu[0][0]);
}
