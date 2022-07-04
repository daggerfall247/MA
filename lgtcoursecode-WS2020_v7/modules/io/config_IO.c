
/*******************************************************************************
 *
 * File config_IO.c
 *
 * Copyright (C) 2013 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines for the reading and writing of the configurations
 *
 * Externally accessible functions:
 *
 * void writeConfig(char *outfile)
 *      Writes the current configuration of gauge fields into the file whose
 *      name is handed to the function via the string `outfile'.
 *
 * void readConfig(char *infile)
 *      Reads the gauge field from the file whose name is handed to the
 *      function via the string infile and stores it into the standard gauge
 *      field.
 *
 *******************************************************************************/

#define CONFIG_IO_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "misc.h"
#include "headers.h"
#include "modules.h"

#define DEBUG_IO 0

#if (MASTER_FIELD == 0)
void writeConfig(char *outfile, int nstart)
{
    FILE *fout = NULL;
    int ii, in, iend;
    long int icheck, icheck2;
    stdint_t lswrite[DIM], info[2];
    double plaq;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("write_config -- in");
    logging("\nWriting configuration %s:\n", outfile);
    t1 = getTime();

#if (DIM == 2)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
#elif (DIM == 3)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
    lswrite[2] = (stdint_t)(LENGS2);
#elif (DIM == 4)
    lswrite[0] = (stdint_t)(LENGT);
    lswrite[1] = (stdint_t)(LENGS1);
    lswrite[2] = (stdint_t)(LENGS2);
    lswrite[3] = (stdint_t)(LENGS3);
#endif

    iend = endianness();
    plaq = plaquette();
    icheck = 0;

    fout = fopen(outfile, "wb");
    error(fout == NULL, "write_config [config_IO.c]",
          "Unable to create output file %s!", outfile);

    info[0] = DIM;
    info[1] = SUN;

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lswrite);
        bswap_double(1, &plaq);
    }

    icheck = fwrite(info, sizeof(stdint_t), 2, fout);
    icheck += fwrite(lswrite, sizeof(stdint_t), DIM, fout);
    icheck += fwrite(&plaq, sizeof(double), 1, fout);
    error(icheck != DIM + 3, "write_config [config_IO.c]", "Write error!");

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "write_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;
    for (in = 0; in < VOL; in++)
    {
        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_sun_dble_array(zw, *pu[in][ii]);

        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }
        icheck += fwrite(buff, sizeof(double), SUNVOL * DIM, fout);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "write_config [config_IO.c]", "Write error!");

    fclose(fout);

    free(buff);

    checkpoint("write_config -- out");
    t2 = getTime();
    logging("\nConfiguration exported (took %.3e sec)!\n", t2 - t1);
}
#elif (MASTER_FIELD == 1)
static long int offset_header;

void writeHeaderToConfig(char *outfile)
{
    FILE *fout = NULL;
    stdint_t lswrite[DIM], info[2];
    double t1, t2;
    int iend, icheck;
    // double plaq;

    checkpoint("write_header_config -- in");
    logging("\nWriting configuration header %s:\n", outfile);
    t1 = getTime();

#if (DIM == 2)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
#elif (DIM == 3)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
    lswrite[2] = (stdint_t)(LENGS2_MF);
#elif (DIM == 4)
    lswrite[0] = (stdint_t)(LENGT_MF);
    lswrite[1] = (stdint_t)(LENGS1_MF);
    lswrite[2] = (stdint_t)(LENGS2_MF);
    lswrite[3] = (stdint_t)(LENGS3_MF);
#endif
    iend = endianness();
    // plaq = plaquette();
    icheck = 0;

    fout = fopen(outfile, "wb");
    error(fout == NULL, "write_config [config_IO.c]", "Unable to create output file %s!", outfile);

    info[0] = DIM;
    info[1] = SUN;

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lswrite);
        // bswap_double(1, &plaq);
    }

    icheck = fwrite(info, sizeof(stdint_t), 2, fout);
    icheck += fwrite(lswrite, sizeof(stdint_t), DIM, fout);
    // icheck += fwrite(&plaq, sizeof(double), 1, fout);

    offset_header = ftell(fout);

    // check 2 statt 3 wegen plaq
    error(icheck != DIM + 2, "write_config [config_IO.c]", "Write error!");
    checkpoint("write_header_config -- out");
    t2 = getTime();
    logging("\nConfiguration exported (took %.3e sec)!\n", t2 - t1);
    fclose(fout);
}

void writeConfig(char *outfile, int nstart)
{
    FILE *fout = NULL;
    int ii, in, iend;
    long int icheck, icheck2;
    long int pos;
    // long int offset_header;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("write_config -- in");
    logging("\nWriting configuration %s at position %d:\n", outfile, nstart);
    t1 = getTime();

    iend = endianness();
    icheck = 0;

    fout = fopen(outfile, "rb+");
    error(fout == NULL, "write_config [config_IO.c]", "Unable to create output file %s!", outfile);

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "write_config [config_IO.c]", "Unable to allocate buffers!");

    // offset_header = 16UL + 4UL * (unsigned long)(DIM);
    for (in = 0; in < VOL; in++)
    {
        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_sun_dble_array(zw, *pu[i[in]][ii]);

        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        pos = index_mf[i[in]];
        pos *= (long int)(sizeof(double) * SUNVOL * DIM);
        fseek(fout, offset_header + pos, SEEK_SET);
        icheck += fwrite(buff, sizeof(double), SUNVOL * DIM, fout);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "write_config [config_IO.c]", "Write error!");

    fclose(fout);

    free(buff);

    checkpoint("write_config -- out");
    t2 = getTime();
    logging("\nConfiguration exported (took %.3e sec)!\n", t2 - t1);
}

void prepareConfig(char *outfile)
{
    int n;
    writeHeaderToConfig(outfile);
    for (n = 0; n < VOL_SL; n++)
    {
        updateArrayOfIndexMF(n);
        writeConfig(outfile, n);
    }
}
#endif

#if (MASTER_FIELD == 0)
void readConfig(char *infile, int nstart)
{
    FILE *fin = NULL;
    int ii, in, iend;
    long int icheck, icheck2;
    int ileng[DIM];
    stdint_t lscheck[DIM], info[2];
    double plaq, plaq0, eps;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("read_config -- in");
    logging("\nReading configuration %s:\n", infile);
    t1 = getTime();

#if (DIM == 2)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
#elif (DIM == 3)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
#elif (DIM == 4)
    ileng[0] = LENGT;
    ileng[1] = LENGS1;
    ileng[2] = LENGS2;
    ileng[3] = LENGS3;
#endif

    error(pu[0][0] == NULL, "read_config [config_IO.c]", "Fields are not allocated!");

    iend = endianness();
    icheck = 0;

    fin = fopen(infile, "rb");
    error(fin == NULL, "read_config [config_IO.c]",
          "Unable to read input file %s!", infile);
    icheck = fread(info, sizeof(stdint_t), 2, fin);
    icheck += fread(lscheck, sizeof(stdint_t), DIM, fin);
    icheck += fread(&plaq0, sizeof(double), 1, fin);
    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lscheck);
        bswap_double(1, &plaq0);
    }
    error(icheck != DIM + 3, "read_config [config_IO.c]", "Read error!");
    error((info[0] != DIM) || (info[1] != SUN), "read_config [config_IO.c]", "Incompatible parameters!");
    for (ii = 0; ii < DIM; ii++)
    {
        error(lscheck[ii] != ileng[ii], "read_config [config_IO.c]", "Incompatible lattice size!");
    }

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "read_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;
    for (in = 0; in < VOL; in++)
    {
        icheck += fread(buff, sizeof(double), SUNVOL * DIM, fin);
        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_dble_array_sun(zw, *pu[in][ii]);
    }

    icheck2 = DIM * SUNVOL * VOL;
    error(icheck != icheck2, "read_config [config_IO.c]", "Read error!");

    fclose(fin);

    plaq = plaquette();
    eps = sqrt((double)(NPLAQ * VOL)) * DBL_EPSILON;
    error(fabs(plaq - plaq0) > eps, "read_config [config_IO.c]", "Plaquette test failed!");

    free(buff);

    checkpoint("read_config -- out");
    t2 = getTime();
    logging("\nConfiguration imported (took %.3e sec)!\n", t2 - t1);
}
#elif (MASTER_FIELD == 1)
void readConfig(char *infile, int nstart)
{
    FILE *fin = NULL;
    int ii, in, iend;
    long int icheck, icheck2;
    long int pos;
    int ileng_mf[DIM];
    stdint_t info[2], lscheck_mf[DIM];
    // double plaq, plaq0, eps;
    double *zw, *buff = NULL;
    double t1, t2;

    checkpoint("read_config -- in");
    logging("\nReading configuration %s at position %d:\n", infile, nstart);
    t1 = getTime();

/*added master-field extents*/
#if (DIM == 2)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
#elif (DIM == 3)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
    ileng_mf[2] = LENGS2_MF;
#elif (DIM == 4)
    ileng_mf[0] = LENGT_MF;
    ileng_mf[1] = LENGS1_MF;
    ileng_mf[2] = LENGS2_MF;
    ileng_mf[3] = LENGS3_MF;
#endif

    error(pu[0][0] == NULL, "read_config [config_IO.c]", "Fields are not allocated!");

    iend = endianness();
    icheck = 0;

    fin = fopen(infile, "rb");
    error(fin == NULL, "read_config [config_IO.c]", "Unable to read input file %s!", infile);
    icheck = fread(info, sizeof(stdint_t), 2, fin);
    icheck += fread(lscheck_mf, sizeof(stdint_t), DIM, fin); /*added master_field extent values*/
    // icheck += fread(&plaq0, sizeof(double), 1, fin);

    if (iend == BIG_ENDIAN)
    {
        bswap_int(2, info);
        bswap_int(DIM, lscheck_mf); /*mf*/
        // bswap_double(1, &plaq0);
    }
    // check 2 statt 3 wegen plaq
    error(icheck != DIM + 2, "read_config [config_IO.c]", "Read error! (Master-field)");
    error((info[0] != DIM) || (info[1] != SUN), "read_config [config_IO.c]", "Incompatible parameters!");
    for (ii = 0; ii < DIM; ii++)
    {
        error(lscheck_mf[ii] != ileng_mf[ii], "read_config [config_IO.c]", "Incompatible master_field size!"); /*added master-field size check*/
    }

    buff = malloc(SUNVOL * DIM * sizeof(double));
    error(buff == NULL, "read_config [config_IO.c]", "Unable to allocate buffers!");

    icheck = 0;

    // offset_header = ftell(fin);

    for (in = 0; in < VOL2; in++)
    {
        pos = index_mf[in];
        pos *= (long int)(SUNVOL * DIM * sizeof(double));
        fseek(fin, offset_header + pos, SEEK_SET); // TODO
        icheck += fread(buff, sizeof(double), SUNVOL * DIM, fin);
        if (iend == BIG_ENDIAN)
        {
            bswap_double(SUNVOL * DIM, buff);
        }

        for (ii = 0, zw = buff; ii < DIM; ii++, zw += SUNVOL)
            mk_dble_array_sun(zw, *pu[in][ii]);
    }

    icheck2 = DIM * SUNVOL * VOL2;
    error(icheck != icheck2, "read_config [config_IO.c]", "Read error!");

    fclose(fin);

    // plaq = plaquette();
    // eps = sqrt((double)(NPLAQ * VOL)) * DBL_EPSILON;
    // error(fabs(plaq - plaq0) > eps, "read_config [config_IO.c]", "Plaquette test failed!");

    free(buff);

    checkpoint("read_config -- out");
    t2 = getTime();
    logging("\nConfiguration imported (took %.3e sec)!\n", t2 - t1);
}

#endif