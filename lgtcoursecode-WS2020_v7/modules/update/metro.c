
/*******************************************************************************
*
* File metro.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to perform Metropolis updates for pure gauge theory.
* 
* Externally accessible functions:
* 
* void staples(int n,int dir,sun_mat *stap)
*      Computes the staples for the link starting at point n in direction
*      dir. The staples matrix is returned via the pointer stap.
*
* double localMetropolisUpdate(int n,int dir,int m)
*        Performs m local Metropolis update steps for the link U_dir(n).
*        On return it hands back the fraction of successful updates.
* 
*******************************************************************************/

#define METRO_C

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"ranlxd.h"
#include"headers.h"
#include"modules.h"


void staples(int n,int dir,sun_mat *stap)
{
   int kk,n1,n2;
   sun_mat un[4];

   sun_zero(un[3]);

   for(kk=0;kk<DIM;kk++)
   {
      if(kk!=dir)
      {
         n1=neib[n][dir];
         n2=neib[n][kk];

         sun_mul_dag(un[1],*pu[n1][kk],*pu[n2][dir]);
         sun_mul_dag(un[0],un[1],*pu[n][kk]);
         sun_add(un[2],un[3],un[0]);

         n2=neib[n][kk+DIM];
	 n1=neib[n2][dir];

         sun_dag(un[0],*pu[n1][kk]);
         sun_mul_dag(un[1],un[0],*pu[n2][dir]);
         sun_mul(un[0],un[1],*pu[n2][kk]);
         sun_add(un[3],un[2],un[0]);
      }
   }

   *stap=un[3];
}


static void proposeNewLink(sun_mat *u)
{
   int i;
   double r[ALGVOL],*a;
   sun_alg X;

   ranlxd(r,ALGVOL);
   a=(double*)(&X);
   for(i=0;i<ALGVOL;i++,a++)
      *a=(1.-2.*r[i]);
   expx(runParams.eps,&X,u);
}

double localMetropolisUpdate(int n,int dir,int m)
{
   int im,iac;
   double lpl,lpl_old,dact,r[1];
   sun_mat stap,upr,zw;

   staples(n,dir,&stap);

   upr=*pu[n][dir];

   for(im=0,iac=0;im<m;im++)
   {
      sun_mul(zw,upr,stap);
      sun_trace(lpl_old,zw);

      proposeNewLink(&upr);
      sun_mul(zw,upr,stap);
      sun_trace(lpl,zw);
      dact=exp(runParams.beta*(lpl-lpl_old)/(double)(SUN));

      ranlxd(r,1);
      if(r[0]<=dact)
      {
         *pu[n][dir]=upr;
         iac++;
      }
      else
         upr=*pu[n][dir];
   }

   return (double)(iac)/m;
}

