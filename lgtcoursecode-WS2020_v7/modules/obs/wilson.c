
/*******************************************************************************
*
* File wilson.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to calculate Wilson loops.
* 
*******************************************************************************/

#define WILSON_C

#include"modules.h"

static int init=0;
static sun_mat *iu1[VOL][DIM],*iu2[VOL][DIM];


static double compute_wloop(int n,int d1,int d2,int l1,int l2)
{
   int ii,nn;
   double tr;
   sun_mat u[2];
   
   sun_unit(u[0]);
   nn=n;
   for(ii=0;ii<l1;ii++)
   {
      u[1]=u[0];
      sun_mul(u[0],u[1],*iu1[nn][d1]);
      nn=neib[nn][d1];
   }
   for(ii=0;ii<l2;ii++)
   {
      u[1]=u[0];
      sun_mul(u[0],u[1],*iu2[nn][d2]);
      nn=neib[nn][d2];
   }
   for(ii=0;ii<l1;ii++)
   {
      u[1]=u[0];
      nn=neib[nn][d1+DIM];
      sun_mul_dag(u[0],u[1],*iu1[nn][d1]);
   }
   for(ii=0;ii<l2;ii++)
   {
      u[1]=u[0];
      nn=neib[nn][d2+DIM];
      sun_mul_dag(u[0],u[1],*iu2[nn][d2]);
   }
   
   sun_trace(tr,u[0]);
   return tr;
}


void measureWilsonLoop(double *w)
{
   int it,ir,nr,nt,nn,inn;
   int in,id;
   int ilin,lbord,nlin,nloop;

   if(init==0)
   {
      allocateGaugeField(iu1);
      allocateGaugeField(iu2);
      init=1;
   }
   nt=((measParams.tf-measParams.ts)/measParams.dt)+1;
   nr=((measParams.rf-measParams.rs)/measParams.dr)+1;
   nloop=nt*nr;
   for(nn=0;nn<nloop;nn++)
      w[nn]=0.;

   copyGaugeField(pu,iu1);
   copyGaugeField(pu,iu2);
   smearing_APE_temporal(measParams.ismt,measParams.smpart,iu1);
   smearing_APE_spatial(measParams.isms,measParams.smpars,iu2);

   nn=0;
   for(it=measParams.ts;it<=measParams.tf;it+=measParams.dt,nn++)
   {
      lbord=latParams.linearExtent[0]-it;
      error(lbord<0,"meas_wilson","Lattice too small to measure loops!");
      nlin=0;
      for(ilin=0;ilin<=lbord;ilin+=it,nlin++)
      {
         inn=nn*nr;
         for(ir=measParams.rs;ir<=measParams.rf;ir+=measParams.dr,inn++)
         {
            for(in=0;in<latParams.volumeOtherDirs[0];in++)
            {
               for (id=1;id<=DIM-1;id++)
               {
                  w[inn]+=compute_wloop(in+ilin*latParams.volumeOtherDirs[0],0,id,it,ir);
               }
            }
         }
      }
      for(ir=0,inn=nn*nr;ir<nr;ir++,inn++)
      {
         w[inn]/=(double)(SUN*latParams.volumeOtherDirs[0]*nlin*(DIM-1));
      }
   }
}
