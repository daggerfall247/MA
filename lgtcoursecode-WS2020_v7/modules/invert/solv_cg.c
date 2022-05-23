
/*******************************************************************************
*
* File solv_cg.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routine to apply the Wilson Dirac operator to a fermion field.
* 
* Externally accessible functions:
* 
* 
* 
*******************************************************************************/

#define SOLV_CG_C

#include<float.h>
#include"headers.h"
#include"modules.h"

static int init=0;
static sun_wferm *z1,*z2,*p,*r;



void allocateFermionFieldsForCG(void)
{
   error(init!=0,"cg [free_cg_fields.c]","CG already initialised!");
   allocateFermionField(&z1);
   allocateFermionField(&z2);
   allocateFermionField(&p);
   allocateFermionField(&r);
   init=1;
}



void deallocateFermionFieldsForCG(void)
{
   error(init==0,"cg [free_cg_fields.c]","CG not initialised!");
   deallocateFermionField(&z1);
   deallocateFermionField(&z2);
   deallocateFermionField(&p);
   deallocateFermionField(&r);
   init=0;
}



int cg(sun_wferm *x,void (*A)(sun_wferm *r,sun_wferm *s),
       void (*Ad)(sun_wferm *r,sun_wferm *s),sun_wferm *b,
       double eps,int nmax)
{
   int n;
   double tol,alp,xi1,xi2;

   error(init==0,"cg [solv_cg.c]","CG not initialised!");

   tol=eps*sqrt(globalSquareNorm(b));

   (*A)(z2,x);
   (*Ad)(z1,z2);
   multiplyByRealAndSum(r,b,-1.,z1);
   xi1=globalSquareNorm(r);
   if(sqrt(xi1)<tol)
      return 0;
   copySpinors(p,r);

   for(n=0;n<nmax;n++)
   {
      (*A)(z2,p);
      (*Ad)(z1,z2);

      xi2=realPartOfScalarProd(p,z1);
      alp=xi1/xi2;
      multiplyByRealAndSum(x,x,alp,p);
      multiplyByRealAndSum(r,r,-alp,z1);

      xi2=globalSquareNorm(x);
      if((100.0*DBL_EPSILON*sqrt(xi2))>tol)
      {
         n=-2;
         break;
      }
      xi2=globalSquareNorm(r);
      if(sqrt(xi2)<tol)
      {
         n++;
         break;
      }

      alp=xi2/xi1;
      xi1=xi2;
      multiplyByRealAndSum(p,r,alp,p);
   }
   if(n==nmax)
      n=-1;

   return n;
}
