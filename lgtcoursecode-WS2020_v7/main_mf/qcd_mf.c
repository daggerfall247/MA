/*******************************************************************************
 *
 * File qcd.c
 *
 * Copyright (C) 2022 Marco Stilger
 * Copyright (C) 2020 Alessandro Sciarra
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Program under construction for lattice gauge theory simulations
 *
 * Syntax: qcd -i [input-filename]
 *
 * For usage instructions see the file ../README.txt
 *
 *******************************************************************************/

#define MAIN_C

#include "ranlxd.h"
#include "modules.h"
#include <stdio.h>
#include <stdbool.h>

int main(int argc, char *argv[])
{
    int n, nw, k;
    bool write;
    int seed, rconf;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12];
    double plaq;

    readInputFile(&seed, out_dir, &rconf, cnfg_file, argc, argv);
    rlxd_init(1, seed);
    initProgram(1);

    setupOutputFiles(runParams.idForOutputFilesName, out_dir);
    printStartupInfo(seed, rconf, cnfg_file);

    initGaugeField(0);

    initGlobalArrays();



    /*
    if(rconf)
       readConfig(cnfg_file,0);

    */
    checkForErrors(1, 0);
    /*
    nw = 0;
    sprintf(cnfg_file, "%s_n%d", CNFG_FILE, nw);
    prepareConfig(cnfg_file);
    */
    
    write = true;
    for (n = 0, nw = 0; n < runParams.numConfs; n++)
    {
        if (write)
        {
            sprintf(cnfg_file, "%s_n%d", CNFG_FILE, nw);
            prepareConfig(cnfg_file);
            write = false;
        }

        plaq = 0.0;
        for (k = 0; k < VOL_SL; k++)
        {
            updateArrayOfIndexMF(k);
            readConfig(cnfg_file, k);
            gaugefieldUpdate(n + 1, 1, SIM_TYPE, SWEEP_TYPE); // adjust number sweeps
            plaq += plaquette();
            writeConfig(cnfg_file, k);
        }
        plaq /= VOL_SL;
        logging("plaq\t%.6f\n",plaq);

        if ((((n - runParams.numThermConfs + 1) % runParams.writeConfsFreq) == 0) && (n > runParams.numThermConfs) && (runParams.writeConfsFreq > 0))
        {
            nw++;
            write = true;
        }

        checkForErrors(1, 0);
    }
    
    releaseGaugeField();

    return 0;
}
