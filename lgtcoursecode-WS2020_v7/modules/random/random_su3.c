
/*******************************************************************************
*
* File random_su3.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to generate a random SU(3) matrix and bosons.
* 
* The routines included are similar to those used by Martin Luescher in
* the DD-HMC code.
* 
* Externally accessible functions:
* 
* void initTwoPi(void)
*      Initialises the local twopi variable.
* 
* int checkTwoPi(void)
*     Checks whether twopi has been initialised.
* 
* void gauss(double *rd,int n)
*      Returns n Gaussian distributed (e^{-r^2}) random numbers in the
*      array rd.
* 
* void su3RandomVector(su3vec *v)
*      Returns a random su3vec v with non-zero norm, including Gaussian
*      distributed random numbers.
*
* void su3RandomMatrix(su3mat *u)
*      Returns a random su3mat u.
*
* void su2RandomMatrix(su2mat *u)
*      Returns a random su2mat u.
*******************************************************************************/

#define RANDOM_SU3_C

#include<float.h>
#include"headers.h"
#include"ranlxd.h"

static int init_twopi=0;
static double twopi;


void initTwoPi(void)
{
   twopi=2.*PI;
   init_twopi=1;
}


int checkTwoPi(void)
{
   return (init_twopi==0);
}


void gauss(double *rd,int n)
{
   int k;
   double ud[2];
   double x1,x2,rho,y1,y2;

   for(k=0;k<n;)
   {
      ranlxd(ud,2);
      x1=ud[0];
      x2=ud[1];
      rho=sqrt(-log(1.0-x1));
      x2*=twopi;
      y1=rho*sin(x2);
      y2=rho*cos(x2);

      rd[k++]=y1;
      if(k<n)
         rd[k++]=y2;
   }
}


void su3RandomVector(su3vec *v)
{
   int i;
   double *r,norm,fact;

   r=(double*)(v);
   norm=0.0;

   while (norm<=DBL_EPSILON)
   {
      gauss(r,6);
      norm=0.0;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=sqrt(norm);
   }

   fact=1.0/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;
}


void su3RandomMatrix(su3mat *u)
{
   int i;
   double *r,norm,fact;
   su3vec *v1,*v2,*v3;

   v1=(su3vec*)(u);
   v2=v1+1;
   v3=v1+2;
   r=(double*)(v3);

   su3RandomVector(v1);
   norm=0.0;

   while (norm<=DBL_EPSILON)
   {
      su3RandomVector(v2);
      su3vec_cross_prod(*v3,*v1,*v2);
      norm=0.0;
      for (i=0;i<6;i++)
         norm+=r[i]*r[i];
      norm=sqrt(norm);
   }

   fact=1.0/norm;
   for (i=0;i<6;i++)
      r[i]*=fact;

   su3vec_cross_prod(*v2,*v3,*v1);
}


void su2RandomMatrix(su2mat *u)
{
   double r[3],zw;

   ranlxd(r,3);
   (*u).c0=1-2.*r[0];
   (*u).c1=1-2.*r[1];
   zw=sqrt(1.-(*u).c1*(*u).c1);
   (*u).c2=zw*cos(twopi*r[2]);
   (*u).c3=zw*sin(twopi*r[2]);
   zw=sqrt(1.-(*u).c0*(*u).c0);
   (*u).c1*=zw;
   (*u).c2*=zw;
   (*u).c3*=zw;
}
