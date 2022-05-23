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
*******************************************************************************/

#define INIT_C

#include"modules.h"
#include<stdlib.h>

#define DEBUG_INIT 0


void initProgram(int itype)
{

   error((runParams.numConfs<=0)||(runParams.decorSteps<0)||(runParams.numThermConfs<0)||
         ((((runParams.numConfs-runParams.numThermConfs)%runParams.writeConfsFreq)!=0)&&(runParams.writeConfsFreq>0)),
         "initProgram [init.c]","Error in run parameters");

   PI=2.0*asin(1.0);

   initTwoPi();
}


void initArrayOfNeighbours(void)
{
   int i,index,dir1,dir2;
   int coord[DIM],newCoord[DIM];
   int volumeFactor[DIM];
   int latticeExtent[DIM],volumeOtherDirs[DIM];

   checkpoint("initArrayOfNeighbours -- in");

#if(DIM==2)
   volumeFactor[0]=LENGS1;
   volumeFactor[1]=1;
   latticeExtent[0]=LENGT;
   latticeExtent[1]=LENGS1;
   volumeOtherDirs[0]=SVOL;
   volumeOtherDirs[1]=LENGT;
#elif(DIM==3)
   volumeFactor[0]=LENGS1*LENGS2;
   volumeFactor[1]=LENGS2;
   volumeFactor[2]=1;
   latticeExtent[0]=LENGT;
   latticeExtent[1]=LENGS1;
   latticeExtent[2]=LENGS2;
   volumeOtherDirs[0]=SVOL;
   volumeOtherDirs[1]=LENGT*LENGS2;
   volumeOtherDirs[2]=LENGT*LENGS1;
#elif(DIM==4)
   volumeFactor[0]=LENGS1*LENGS2*LENGS3;
   volumeFactor[1]=LENGS2*LENGS3;
   volumeFactor[2]=LENGS3;
   volumeFactor[3]=1;
   latticeExtent[0]=LENGT;
   latticeExtent[1]=LENGS1;
   latticeExtent[2]=LENGS2;
   latticeExtent[3]=LENGS3;
   volumeOtherDirs[0]=SVOL;
   volumeOtherDirs[1]=LENGT*LENGS2*LENGS3;
   volumeOtherDirs[2]=LENGT*LENGS1*LENGS3;
   volumeOtherDirs[3]=LENGT*LENGS1*LENGS2;
#endif

   for(dir1=0;dir1<DIM;dir1++)
   {
      latParams.linearExtent[dir1]=latticeExtent[dir1];
      latParams.volumeOtherDirs[dir1]=volumeOtherDirs[dir1];
      latParams.volumeFactor[dir1]=volumeFactor[dir1];
   }

   for(i=0;i<VOL;i++)
   {
      // Compute coordinates in all directions for site index i
      for(dir1=0,index=i;dir1<DIM;dir1++)
      {
         coord[dir1]=index/volumeFactor[dir1];
         index=index%volumeFactor[dir1];
      }
      for(dir1=0;dir1<DIM;dir1++)
      {
         for(dir2=0;dir2<DIM;dir2++)
            newCoord[dir2]=coord[dir2];

         // Compute coordinates of nearest neighbour in positive directions
         newCoord[dir1]=(coord[dir1]+1)%latticeExtent[dir1];
         for(dir2=0,neib[i][dir1]=0;dir2<DIM;dir2++)
            neib[i][dir1]+=newCoord[dir2]*volumeFactor[dir2];

         // Compute coordinates of nearest neighbour in negative directions
         newCoord[dir1]=(coord[dir1]+latticeExtent[dir1]-1)%latticeExtent[dir1];
         for(dir2=0,neib[i][dir1+DIM]=0;dir2<DIM;dir2++)
            neib[i][dir1+DIM]+=newCoord[dir2]*volumeFactor[dir2];
      }
   }

#if(DEBUG_INIT==1)
   logging("\nix");
   for(dir1=0;dir1<DIM;dir1++)
      logging("\tn%d\tn%d",dir1,dir1+DIM);
   for(i=0;i<VOL;i++)
   {
      logging("\n%d",i);
      for(dir1=0;dir1<DIM;dir1++)
         logging("\t%d\t%d",neib[i][dir1],neib[i][dir1+DIM]);
   }
   logging("\n");
#endif

   checkpoint("initArrayOfNeighbours -- out");
   return;
}


void initGaugeField(int flag)
{
   sun_mat *iu=NULL,*zu,ini;
   int ind,dir;

   checkpoint("init_gauge");

   error((flag!=0)&&(flag!=1),"init_gauge [start.c]",
         "Wrong flag! Cannot initiate gauge field!");

   iu=malloc(DIM*VOL*sizeof(sun_mat));
   error(iu==NULL,"init_gauge [start.c]",
         "Unable to allocate gauge field array!");

   if(flag==0)
   {
      sun_unit(ini);
   }
   for(ind=0,zu=iu;ind<VOL;ind++)
   {
      for(dir=0;dir<DIM;dir++,zu++)
      {
         pu[ind][dir]=zu;
         if(flag==0)
         {
            *pu[ind][dir]=ini;
         }
         else if(flag==1)
         {
#if(SUN==3)
            su3RandomMatrix(pu[ind][dir]);
#elif(SUN==2)
            su2RandomMatrix(pu[ind][dir]);
#endif
         }
      }
   }
}


void releaseGaugeField(void)
{
   free(pu[0][0]);
}


void allocateGaugeField(sun_mat *u[VOL][DIM])
{
   sun_mat *iu=NULL,*zu;
   int ind,dir;

   checkpoint("allocate_gauge");

   iu=malloc(DIM*VOL*sizeof(sun_mat));
   error(iu==NULL,"allocate_gauge [start.c]",
         "Unable to allocate gauge field array!");

   for(ind=0,zu=iu;ind<VOL;ind++)
   {
      for(dir=0;dir<DIM;dir++,zu++)
      {
         u[ind][dir]=zu;
      }
   }
}


void deallocateGaugeField(sun_mat *u[VOL][DIM])
{
   free(u[0][0]);
}


void copyGaugeField(sun_mat *u1[VOL][DIM],sun_mat *u2[VOL][DIM])
{
   int ii,jj;
   for(ii=0;ii<VOL;ii++)
      for(jj=0;jj<DIM;jj++)
         *u2[ii][jj]=*u1[ii][jj];
}


void allocateFermionField(sun_wferm **f)
{
   sun_wferm *s;
   s=malloc(VOL*sizeof(sun_wferm));
   error(s==NULL,"allocateFermionField [init.c]",
        "Unable to allocate fermion field array!");
   *f=s;
}


void deallocateFermionField(sun_wferm **f)
{
    free(*f);
}
