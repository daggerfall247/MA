
/*******************************************************************************
*
* File spin_alg.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routine needed for spinor algebra.
* 
* Externally accessible functions:
* 
* void setAllSpinorsToZero(sun_wferm *f)
*      Sets the all entries of the spinor pointed to by f to 0.
* 
* void copySpinors(sun_wferm *f2,sun_wferm *f1)
*      Copies the spinor f1 into f2.
* 
* void multiplyByRealAndSum(sun_wferm *f3,sun_wferm *f2,double a,sun_wferm *f1)
*      Computes f3=f2+a*f1 where f1,f2,f3 are the spinors pointed to by
*      the associated pointers and a is a real number.
* 
* double globalSquareNorm(sun_wferm *f)
*     Returns |f|^2, normalised to 1 for a unit spinor (i.e. the 'volume'
*     factor DIM*SUN*VOL is divided out).
* 
* double globalSum(sun_wferm *f)
*     Returns the sum of all components, (real and imaginary ones)
*     normalised to 1 for a unit spinor.
* 
* double realPartOfScalarProd(sun_wferm *f2,sun_wferm *f1)
*         Returns the real part of the scalar product (f2,f1) normalised
*         to 1 when f2 and f1 are unit spinors.
*
*******************************************************************************/

#define SPIN_ALG_C

#include"modules.h"



void setAllSpinorsToZero(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+2*DIM*SUN*VOL;
    for(;s<sf;s+=1)
       *s=0.;
}


void copySpinors(sun_wferm *f2,sun_wferm *f1)
{
   double *s1,*s2,*sf;
   s1=(double*)(f1);
   s2=(double*)(f2);
   sf=s1+2*DIM*SUN*VOL;
   for(;s1<sf;s1+=1,s2+=1)
      *s2=*s1;
}


void multiplyByRealAndSum(sun_wferm *f3,sun_wferm *f2,double a,sun_wferm *f1)
{
   double *s1,*s2,*s3,*sf;
   s1=(double*)(f1);
   s2=(double*)(f2);
   s3=(double*)(f3);
   sf=s1+2*DIM*SUN*VOL;
   for(;s1<sf;s1+=1,s2+=1,s3+=1)
      *s3=*s2+a*(*s1);
}


double globalSquareNorm(sun_wferm *f)
{
   double norm;
   complex *s,*sf;
   s=(complex*)(f);
   sf=s+DIM*SUN*VOL;
   for(norm=0.;s<sf;s+=1)
      norm+=(*s).re*(*s).re+(*s).im*(*s).im;
   return norm/(double)(DIM*SUN*VOL);
}


double globalSum(sun_wferm *f)
{
   double sum;
   complex *s,*sf;
   s=(complex*)(f);
   sf=s+DIM*SUN*VOL;
   for(sum=0.;s<sf;s+=1)
   {
      sum+=(*s).re;
      sum+=(*s).im;
   }
   return sum/(double)(DIM*SUN*VOL);
}


complex scalarProd(sun_wferm *f2,sun_wferm *f1)
{
   complex z,prod;
   complex *s1,*s2,*sf;

   s1=(complex*)(f1);
   s2=(complex*)(f2);
   sf=s1+4*SUN*VOL;
   prod.re=0.;
   prod.im=0.;
   for(;s1<sf;s1+=1,s2+=1)
   {
      compl_mult_sn(z,*s2,*s1);
      compl_selfadd(prod,z);
   }
   compl_realdiv_single(prod,(double)(4*SUN*VOL));

   return prod;
}


double realPartOfScalarProd(sun_wferm *f2,sun_wferm *f1)
{
   double prod,z;
   complex *s1,*s2,*sf;

   s1=(complex*)(f1);
   s2=(complex*)(f2);
   sf=s1+4*SUN*VOL;
   for(prod=0.;s1<sf;s1+=1,s2+=1)
   {
      compl_mult_sn_re(z,*s2,*s1);
      prod+=z;
   }

   return prod/(double)(DIM*SUN*VOL);
}
