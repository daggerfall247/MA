
/*******************************************************************************
*
* File qcd.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program for
*
* Syntax: meas -i [input-filename]
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"ranlxd.h"
#include"modules.h"


void setGaugeFieldToCold(void)
{
   int n,dir;
   for(n=0;n<VOL;n++)
      for(dir=0;dir<DIM;dir++)
         sun_unit(*pu[n][dir]);
}


int main(int argc,char *argv[])
{
   int n,ir,it,nr,nt,il;
   int seed;
   int sconf,fconf,econf;
   char out_dir[NAME_SIZE];
   char cnfg_file[FULL_PATH_SIZE+12];
   double *wil;

   readInputFileForMeas(&seed,&sconf,&fconf,&econf,out_dir,argc,argv);
   rlxd_init(1,seed);

   setupOutputFilesForMeas(runParams.idForOutputFilesName,out_dir);
   initArrayOfNeighbours();
   initGaugeField(1);

   nt=((measParams.tf-measParams.ts)/measParams.dt)+1;
   nr=((measParams.rf-measParams.rs)/measParams.dr)+1;
   if(runParams.mwil)
      wil=malloc(nr*nt*sizeof(double));
   if(runParams.mcorrs)
      allocateFermionFieldsFor2ptFunctions(measParams.stype);

   for(n=sconf;n<=fconf;n+=econf)
   {
      sprintf(cnfg_file,"%s_n%d",CNFG_FILE,n);
      logging("\nMeasuring on configuration:\n %s\n",cnfg_file);
      readConfig(cnfg_file);

      if(runParams.mwil)
      {
         measureWilsonLoop(wil);
         for(it=0,il=0;it<nt;it++)
         {
            for(ir=0;ir<nr;ir++,il++)
            {
               logging("WIL %d %d\t%.8e\n",
                  measParams.ts+it*measParams.dt,measParams.rs+ir*measParams.dr,
                  wil[il]);
            }
         }
      }
      if(runParams.mcorrs)
      {
         measure2ptFunctions(measParams.source,measParams.stype,measParams.n2pt,measParams.g1,measParams.g2);
      }
   }
   if(runParams.mwil)
      free(wil);
   if(runParams.mcorrs)
   {
      deallocateFermionFieldsFor2ptFunctions();
      free(measParams.g1);
      free(measParams.g2);
   }


   releaseGaugeField();

   return 0;
}
