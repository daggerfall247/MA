
/*******************************************************************************
*
* File test.c

* Copyright (C) 2021 Alessandro Sciarra
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program under construction for testing lattice gauge theory code
*
* Syntax: test -i <input-filename>
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include "ranlxd.h"
#include "modules.h"
#include "test_utils.h"
#include <assert.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int seed, rconf;
    char out_dir[NAME_SIZE];
    char cnfg_file[FULL_PATH_SIZE + 12];
    readInputFile(&seed, out_dir, &rconf, cnfg_file, argc, argv);
    setupOutputFiles(runParams.idForOutputFilesName,out_dir);
    initArrayOfNeighbours();
    initGaugeField(0);
    printf("%f\n", plaquette());
    printf("test1\n");
    readConfig("test.bin", 0);
    printf("%f\n", plaquette());

    printf("test2\n");
    //printf("%f %f %f %f", (*pu[0][0]).c0, (*pu[0][0]).c1, (*pu[0][0]).c2, (*pu[0][0]).c3);
    releaseGaugeField();
    /*runAllTests();*/

    return 0;
}
