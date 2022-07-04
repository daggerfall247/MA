
/*******************************************************************************
 *
 * File inp_IO.c
 *
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Includes the routines for the reading of the infile and the setup of the
 * files, as well as some startup writeouts.
 *
 * Externally accessile functions:
 *
 * void readInputFile(int *iseed,char *dir,int *rconf,
 *               char *cfile,int argc,char *argv[])
 *      Reads the input file and stores the associated parameters in the
 *      struct 'runParameters runp' (defined in headers.h). Some infos are
 *      also pased back to the calling routine (typically the main) via
 *      the arguments of the function.
 *
 *      For parameters in runParameters see the documentation.
 *      iseed : seed for random number generator
 *      dir: directory for output
 *      rconf: Start from a given configuration?
 *      cfile: Filename of that configuration
 *
 * void setupOutputFiles(int id,char *dir)
 *      Setup of the main files for output.
 *
 *      id: run id
 *      dir: output diectory
 *
 * void printStartupInfo(int iseed,int rconf,char *cfile)
 *      Printing the startup information into the logfile. for identification
 *      of run parameters etc..
 *
 *      iseed: seed of randum number generator
 *      rconf,cfile: as above
 *
 *******************************************************************************/

#define INP_IO_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include "modules.h"

void readInputFile(int *iseed, char *dir, int *rconf, char *cfile, int argc, char *argv[])
{
    FILE *inf = NULL, *ftest = NULL;
    int i, ii, ifail = 0;
    int id, nconf, nstep, ntherm, wcnfg;
    double beta, eps;

    ftest = fopen("SETUP_ERRORS", "r");
    if (ftest != NULL)
    {
        fclose(ftest);
        remove("SETUP_ERRORS");
    }
    sprintf(LOG_FILE, "SETUP_ERRORS");

    for (i = 1, ii = 0; i < argc; i++)
        if (strcmp(argv[i], "-c") == 0)
            ii = i + 1;
    if (ii != 0)
    {
        strcpy(cfile, argv[ii]);
        inf = fopen(cfile, "r");
        error(inf == NULL, 0, "readInputFile [inp_IO.c]",
              "Incorrect config file!\n"
              "Syntax: IND-G_pg -i <input file> [-c <input config>]");
        fclose(inf);
        inf = NULL;
        *rconf = 1;
    }
    else
        *rconf = 0;

    for (i = 1, ii = 0; i < argc; i++)
        if (strcmp(argv[i], "-i") == 0)
            ii = i + 1;
    error(ii == 0, "readInputFile [inp_IO.c]",
          "Syntax: IND-G_pg -i <input file> [-c <input config>]");
    inf = fopen(argv[ii], "r");
    error(inf == NULL, "readInputFile [inp_IO.c]",
          "Unable to open input file!");

    ifail += fscanf(inf, "outFileId     %d\n", &id);
    ifail += fscanf(inf, "outputDir     %s\n", dir);
    ifail += fscanf(inf, "beta          %lf\n", &beta);
    ifail += fscanf(inf, "epsilon       %lf\n", &eps);
    ifail += fscanf(inf, "nconf         %d\n", &nconf);
    ifail += fscanf(inf, "ntherm        %d\n", &ntherm);
    ifail += fscanf(inf, "ndecorr       %d\n", &nstep);
    ifail += fscanf(inf, "writeConfFreq %d\n", &wcnfg);
    ifail += fscanf(inf, "seed          %d\n", iseed);
    error(ifail != 9, "readInputFile [inp_IO.c]",
          "Unable to read some fields in input file");

    fclose(inf);

    runParams.idForOutputFilesName = id;
    runParams.beta = beta;
    runParams.eps = eps;
    runParams.numConfs = nconf;
    runParams.decorSteps = nstep;
    runParams.numThermConfs = ntherm;
    runParams.writeConfsFreq = wcnfg;
}

void setupOutputFiles(int id, char *dir)
{
    FILE *fset = NULL;
    char base[NAME_SIZE], check[FULL_PATH_SIZE];

#if (MASTER_FIELD == 1)
#if (DIM == 2)
    sprintf(base, "%dx%d_SU%d_b%.6f_id%d", LENGT_MF, LENGS1_MF, SUN, runParams.beta, id);
#elif (DIM == 3)
    sprintf(base, "%dx%dx%d_SU%d_b%.6f_id%d", LENGT_MF, LENGS1_MF, LENGS2_MF, SUN, runParams.beta, id);
#elif (DIM == 4)
    sprintf(base, "%dx%dx%dx%d_SU%d_b%.6f_id%d", LENGT_MF, LENGS1_MF, LENGS2_MF, LENGS3_MF, SUN, runParams.beta, id);
#endif
#else
#if (DIM == 2)
    sprintf(base, "%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1, SUN, runParams.beta, id);
#elif (DIM == 3)
    sprintf(base, "%dx%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1, LENGS2, SUN, runParams.beta, id);
#elif (DIM == 4)
    sprintf(base, "%dx%dx%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1, LENGS2, LENGS3, SUN, runParams.beta, id);
#endif
#endif

    sprintf(check, "%s/%s.log", dir, base);
    fset = fopen(check, "r");
    error(fset != NULL, "setupOutputFiles [inp_IO.c]",
          "Try to overwrite existing .log file!");
    sprintf(check, "%s/%s.out", dir, base);
    fset = fopen(check, "r");
    error(fset != NULL, "setupOutputFiles [inp_IO.c]",
          "Try to overwrite existing .out file!");

    fset = fopen(LOG_FILE, "r");
    if (fset != NULL)
    {
        printf("SETUP_ERRORS file was not empty!!!\n");
        printf("Aborting!\n");
        exit(-1);
    }

    sprintf(LOG_FILE, "%s/%s.log", dir, base);
    sprintf(OUT_FILE, "%s/%s.out", dir, base);
    sprintf(CNFG_FILE, "%s/%s", dir, base);

    fset = fopen(LOG_FILE, "w");
    error(fset == NULL, "setupOutputFiles [inp_IO.c]",
          "Unable to create .log file!");
    fclose(fset);
    fset = fopen(OUT_FILE, "w");
    error(fset == NULL, "setupOutputFiles [inp_IO.c]",
          "Unable to create .out file!");
    fclose(fset);
}

void printStartupInfo(int iseed, int rconf, char *cfile)
{
    FILE *flog = NULL;

    flog = fopen(LOG_FILE, "ab");
    error(flog == NULL, "printStartupInfo [inp_IO.c]",
          "Unable to open .log file!");

    fprintf(flog, "*******************************\n");
    fprintf(flog, "LATTICE QCD CODE\n");
    fprintf(flog, "*******************************\n");
    fprintf(flog, "Global parameters:\n");
    fprintf(flog, "*******************************\n");
#if(MASTER_FIELD==1)
    fprintf(flog, "Master_field simulation:\n");
#if (DIM == 2)
    fprintf(flog, "Using a two-dimensional master-field:\n");
    fprintf(flog, "master-field  :\t%dx%d\n", LENGT_MF, LENGS1_MF);
    fprintf(flog, "with %d sublattices\n", FRACT*FRACS1);
    fprintf(flog, "of size       :\t%dx%d\n", LENGT, LENGS1);
#elif (DIM == 3)
    fprintf(flog, "Using a three-dimensional master-field:\n");
    fprintf(flog, "master-field  :\t%dx%dx%d\n", LENGT_MF, LENGS1_MF, LENGS2_MF);
    fprintf(flog, "with %d sublattices\n", FRACT*FRACS1*FRACS2);
    fprintf(flog, "of size       :\t%dx%dx%d\n", LENGT, LENGS1, LENGS2);
#else
    fprintf(flog, "Using a four -dimensional master-field:\n");
    fprintf(flog, "master-field  :\t%dx%dx%dx%d\n", LENGT_MF, LENGS1_MF, LENGS2_MF, LENGS3_MF);
    fprintf(flog, "with %d sublattices\n", FRACT*FRACS1*FRACS2*FRACS3);
    fprintf(flog, "of size       :\t%dx%dx%dx%d\n", LENGT, LENGS1, LENGS2, LENGS3);
#endif
#else
#if (DIM == 2)
    fprintf(flog, "Using a two-dimensional lattice:\n");
    fprintf(flog, "latt  :\t%dx%d\n", LENGT, LENGS1);
#elif (DIM == 3)
    fprintf(flog, "Using a three-dimensional lattice:\n");
    fprintf(flog, "latt  :\t%dx%dx%d\n", LENGT, LENGS1, LENGS2);
#else
    fprintf(flog, "Using a four -dimensional lattice:\n");
    fprintf(flog, "latt  :\t%dx%dx%dx%d\n", LENGT, LENGS1, LENGS2, LENGS3);
#endif
#endif
    fprintf(flog, "beta :\t%.6f\n", runParams.beta);
    fprintf(flog, "epsilon :\t%.6f\n", runParams.eps);
    fprintf(flog, "nconf :\t%d\n", runParams.numConfs);
    fprintf(flog, "ndecorr :\t%d\n", runParams.decorSteps);
    fprintf(flog, "ntherm :\t%d\n", runParams.numThermConfs);
    fprintf(flog, "writeConfFreq :\t%d\n", runParams.writeConfsFreq);
    fprintf(flog, "seed :\t%d\n", iseed);
    fprintf(flog, "*******************************\n");
    fprintf(flog, "Inline measurements:\n");
    fprintf(flog, "None at the moment!\n");
    fprintf(flog, "*******************************\n");
    fprintf(flog, "Start from \"%s\" configuration:\n", cfile);
    fprintf(flog, "*******************************\n");

    fclose(flog);
}

void readInputFileForMeas(int *iseed, int *sconf, int *fconf, int *econf, char *dir, int argc, char *argv[])
{
    FILE *inf = NULL, *ftest = NULL;
    int i, ii, ifail = 0;
    double kappa;

    ftest = fopen("SETUP_ERRORS", "r");
    if (ftest != NULL)
    {
        fclose(ftest);
        remove("SETUP_ERRORS");
    }
    sprintf(LOG_FILE, "SETUP_ERRORS");

    for (i = 1, ii = 0; i < argc; i++)
        if (strcmp(argv[i], "-i") == 0)
            ii = i + 1;
    error(ii == 0, "readInputFileForMeas [inp_IO.c]",
          "Syntax: meas -i <input file> [-c <input config>]");
    inf = fopen(argv[ii], "r");
    error(inf == NULL, "readInputFileForMeas [inp_IO.c]",
          "Unable to open input file!");

    ifail += fscanf(inf, "outFileId      %d\n", &runParams.idForOutputFilesName);
    ifail += fscanf(inf, "odir           %s\n", dir);
    ifail += fscanf(inf, "beta           %lf\n", &runParams.beta);
    ifail += fscanf(inf, "seed           %d\n", iseed);
    ifail += fscanf(inf, "sconf          %d\n", sconf);
    ifail += fscanf(inf, "fconf          %d\n", fconf);
    ifail += fscanf(inf, "econf          %d\n", econf);
    ifail += fscanf(inf, "mwil           %d\n", &runParams.mwil);
    ifail += fscanf(inf, "tLoopExtent    %d %d %d\n", &measParams.ts, &measParams.tf, &measParams.dt);
    ifail += fscanf(inf, "tSmear         %d %lf\n", &measParams.ismt, &measParams.smpart);
    ifail += fscanf(inf, "rLoopExtent    %d %d %d\n", &measParams.rs, &measParams.rf, &measParams.dr);
    ifail += fscanf(inf, "rSmear         %d %lf\n", &measParams.isms, &measParams.smpars);
    ifail += fscanf(inf, "mcorrs         %d\n", &runParams.mcorrs);
    ifail += fscanf(inf, "kappa          %lf\n", &kappa);
    ifail += fscanf(inf, "source         %d %d %d %d\n", measParams.source, measParams.source + 1, measParams.source + 2, measParams.source + 3);
    ifail += fscanf(inf, "stype          %d\n", &measParams.stype);
    ifail += fscanf(inf, "nmax           %d\n", &measParams.nmax);
    ifail += fscanf(inf, "prec           %lf\n", &measParams.eps);
    ifail += fscanf(inf, "n2pt           %d\n", &measParams.n2pt);
    error(ifail != 28, "readInputFileForMeas [inp_IO.c]",
          "Unable to read some fields in input file");
    measParams.g1 = malloc(measParams.n2pt * sizeof(int));
    measParams.g2 = malloc(measParams.n2pt * sizeof(int));
    ifail = 0;
    for (ii = 0; ii < measParams.n2pt; ii++)
        ifail += fscanf(inf, "gammas      %d %d\n", measParams.g1 + ii, measParams.g2 + ii);
    error(ifail != 2 * measParams.n2pt, "readInputFileForMeas [inp_IO.c]",
          "Unable to read some gammas fields in input file");

    fclose(inf);

    runParams.mass = 1. / (2 * kappa) - 4;
}

void setupOutputFilesForMeas(int id, char *dir)
{
    FILE *fset = NULL;
    char base[NAME_SIZE];

    if (DIM == 2)
        sprintf(base, "%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1,
                SUN, runParams.beta, id);
    else if (DIM == 3)
        sprintf(base, "%dx%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1, LENGS2,
                SUN, runParams.beta, id);
    else
        sprintf(base, "%dx%dx%dx%d_SU%d_b%.6f_id%d", LENGT, LENGS1, LENGS2, LENGS3,
                SUN, runParams.beta, id);

    fset = fopen("SETUP_ERRORS", "r");
    if (fset != NULL)
    {
        fclose(fset);
        remove("SETUP_ERRORS");
    }
    sprintf(LOG_FILE, "SETUP_ERRORS");

    DIR *openDir = opendir(dir);
    error(openDir == NULL, "setupOutputFilesForMeas [inp_IO.c]",
          "Unable to open %s directory!, dir");
    if (openDir != NULL)
        closedir(openDir); /* Directory exists. */

    sprintf(LOG_FILE, "%s/%s.meas", dir, base);
    sprintf(CNFG_FILE, "%s/%s", dir, base);

    fset = fopen(LOG_FILE, "w");
    error(fset == NULL, "setupOutputFilesForMeas [inp_IO.c]",
          "Unable to create .meas file!");
    fclose(fset);
}
