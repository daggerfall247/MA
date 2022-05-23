
/*******************************************************************************
*
* File smearing.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to perform APE smearing of a gauge configuration.
* 
* Externally accessible functions:
* 
* void smearing_APE_spatial(int ns,double sp,sun_mat *u[VOL][DIM])
*      Performs ns smearing levels of APE smearing for the spatial
*      links included in the gauge field u. Note, in the course
*      of the smearing no temporal staples are used. The smearing
*      parameter is handed over to the routine in sp.
* 
* void smearing_APE_temporal(int ns,double sp,sun_mat *u[VOL][DIM])
*      Performs ns smearing levels of APE smearing for the temporal
*      links included in the gauge field u. Note, in the course
*      of the smearing no spatial links are changed. The smearing
*      parameter is handed over to the routine in sp.
* 
*******************************************************************************/

#define SMEARING_C

#include<stdio.h>
#include"headers.h"
#include"modules.h"

static int init=0;
static sun_mat *zu[VOL][DIM];


static void staples_smearing(int n,int dir,int dt,sun_mat *stap)
{
   int kk,n1,n2;
   sun_mat un[4];

   sun_zero(un[3]);

   for(kk=0;kk<DIM;kk++)
   {
      if((kk!=dir)&&(kk!=dt))
      {
         n1=neib[n][dir];
         n2=neib[n][kk];

         sun_mul(un[1],*zu[n][kk],*zu[n2][dir]);
         sun_mul_dag(un[0],un[1],*zu[n1][kk]);
         sun_add(un[2],un[3],un[0]);

         n2=neib[n][kk+DIM];
         n1=neib[n2][dir];

         sun_dag(un[0],*zu[n2][kk]);
         sun_mul(un[1],un[0],*zu[n2][dir]);
         sun_mul(un[0],un[1],*zu[n1][kk]);
         sun_add(un[3],un[2],un[0]);
      }
   }

   *stap=un[3];
}



void smearing_APE_spatial(int ns,double sp,sun_mat *u[VOL][DIM])
{
   int ii,jj,n;
   sun_mat stap;

   if(init==0)
   {
      allocateGaugeField(zu);
      init=1;
   }

   for(n=0;n<ns;n++)
   {
      copyGaugeField(u,zu);
      for(ii=0;ii<VOL;ii++)
      {
         for(jj=1;jj<DIM;jj++)
         {
            staples_smearing(ii,jj,0,&stap);
            sun_dble_mul(stap,sp/(2*(DIM-2)));
            sun_dble_mul(*u[ii][jj],1.-sp);
            sun_add(*u[ii][jj],*u[ii][jj],stap);
#if(SUN==2)
            project_to_su2(u[ii][jj]);
#elif(SUN==3)
            project_to_su3(u[ii][jj]);
#endif
         }
      }
   }
}


void smearing_APE_temporal(int ns,double sp,sun_mat *u[VOL][DIM])
{
   int ii,n;
   sun_mat stap;

   if(init==0)
   {
      allocateGaugeField(zu);
      init=1;
   }
   
   for(n=0;n<ns;n++)
   {
      copyGaugeField(u,zu);
      for(ii=0;ii<VOL;ii++)
      {
         staples_smearing(ii,0,0,&stap);
         sun_dble_mul(stap,sp/(2*(DIM-1)));
         sun_dble_mul(*u[ii][0],1.-sp);
         sun_add(*u[ii][0],*u[ii][0],stap);
#if(SUN==2)
         project_to_su2(u[ii][0]);
#elif(SUN==3)
         project_to_su3(u[ii][0]);
#endif
      }
   }
}
