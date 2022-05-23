
/*******************************************************************************
*
* File exp_fct.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to determine the exponential function exp(X),
* where X is an element of the Lie algebra of SU(N), correctly to O(X^2).
* 
* The routines included are similar to those used by Martin Luescher in
* the DD-HMC code.
* 
* Externally accessible functions:
* 
* void expx(double eps,su3_alg_dble *X,su3_dble *u)
*      Replaces u by exp(eps*X)*u, where "exp" is the approximate SU(N)
*      exponential function.
* 
* Note: The implementation in this version should only be used for the
*       leapfrog integrator in the HMC for other integrator schemes it
*       might not be accurate enough!
*       -> In that case: Use the Cayley-Hamilton exponential function!
*
*******************************************************************************/

#define EXP_FCT_C

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include"headers.h"
#include"modules.h"



#if(SUN==3)

static double s1[4],s2[4],s3[4];
static double r1;
static complex z1,z2;



static void su2_rotate(double *s,su3vec *v1,su3vec *v2)
{
   z1.re=
      s[0]*(*v1).c1.re-s[1]*(*v2).c1.im+s[2]*(*v2).c1.re-s[3]*(*v1).c1.im;
   z1.im=
      s[0]*(*v1).c1.im+s[1]*(*v2).c1.re+s[2]*(*v2).c1.im+s[3]*(*v1).c1.re;
   z2.re=
      s[0]*(*v2).c1.re-s[1]*(*v1).c1.im-s[2]*(*v1).c1.re+s[3]*(*v2).c1.im;
   z2.im=
      s[0]*(*v2).c1.im+s[1]*(*v1).c1.re-s[2]*(*v1).c1.im-s[3]*(*v2).c1.re;
   (*v1).c1=z1;
   (*v2).c1=z2;

   z1.re=
      s[0]*(*v1).c2.re-s[1]*(*v2).c2.im+s[2]*(*v2).c2.re-s[3]*(*v1).c2.im;
   z1.im=
      s[0]*(*v1).c2.im+s[1]*(*v2).c2.re+s[2]*(*v2).c2.im+s[3]*(*v1).c2.re;
   z2.re=
      s[0]*(*v2).c2.re-s[1]*(*v1).c2.im-s[2]*(*v1).c2.re+s[3]*(*v2).c2.im;
   z2.im=
      s[0]*(*v2).c2.im+s[1]*(*v1).c2.re-s[2]*(*v1).c2.im-s[3]*(*v2).c2.re;
   (*v1).c2=z1;
   (*v2).c2=z2;

   z1.re=
      s[0]*(*v1).c3.re-s[1]*(*v2).c3.im+s[2]*(*v2).c3.re-s[3]*(*v1).c3.im;
   z1.im=
      s[0]*(*v1).c3.im+s[1]*(*v2).c3.re+s[2]*(*v2).c3.im+s[3]*(*v1).c3.re;
   z2.re=
      s[0]*(*v2).c3.re-s[1]*(*v1).c3.im-s[2]*(*v1).c3.re+s[3]*(*v2).c3.im;
   z2.im=
      s[0]*(*v2).c3.im+s[1]*(*v1).c3.re-s[2]*(*v1).c3.im-s[3]*(*v2).c3.re;
   (*v1).c3=z1;
   (*v2).c3=z2;   
}



void expx(double t,su3alg *X,su3mat *u)
{
   su3vec *v1,*v2,*v3;

   s1[1]=t*(*X).c4;
   s1[2]=t*(*X).c3;
   s1[3]=t*(*X).c1;

   r1=s1[1]*s1[1]+s1[2]*s1[2]+s1[3]*s1[3];   
   r1=1.0/(2.0+0.125*r1);
   s1[0]=4.0*r1-1.0;
   s1[1]*=r1;
   s1[2]*=r1;
   s1[3]*=r1;   
   
   s2[1]=t*(*X).c6;
   s2[2]=t*(*X).c5;
   s2[3]=t*(*X).c2;   

   r1=s2[1]*s2[1]+s2[2]*s2[2]+s2[3]*s2[3];   
   r1=1.0/(2.0+0.125*r1);
   s2[0]=4.0*r1-1.0;
   s2[1]*=r1;
   s2[2]*=r1;
   s2[3]*=r1;  
   
   s3[1]=t*(*X).c8;
   s3[2]=t*(*X).c7;
   s3[3]=t*((*X).c2-(*X).c1);
   
   r1=s3[1]*s3[1]+s3[2]*s3[2]+s3[3]*s3[3];   
   r1=1.0/(1.0+0.25*r1);
   s3[0]=2.0*r1-1.0;
   s3[1]*=r1;
   s3[2]*=r1;
   s3[3]*=r1;

   v1=(su3vec*)(u);
   v2=v1+1;
   v3=v1+2;
   
   su2_rotate(s1,v1,v2);
   su2_rotate(s2,v1,v3);   
   su2_rotate(s3,v2,v3);
   su2_rotate(s2,v1,v3);
   su2_rotate(s1,v1,v2);
}

#elif(SUN==2)

void expx(double t,su2alg *X,su2mat *u)
{
   double zw;
   su2mat etx,uz;
   
   zw=t*t*((*X).c1*(*X).c1+(*X).c2*(*X).c2+(*X).c3*(*X).c3)/4.;
   etx.c0=1.-zw;
   zw=1./(1.+zw);
   etx.c0*=zw;
   etx.c1=zw*t*(*X).c1;
   etx.c2=zw*t*(*X).c2;
   etx.c3=zw*t*(*X).c3;
   su2_mat_mul(uz,etx,*u);
   *u=uz;
}

#endif
