
/*******************************************************************************
*
* File dirac_wil.c
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
* applyWilsonDiracOperator(sun_wferm *r,sun_wferm *s)
* 
* apply_dirac_wil_conj(sun_wferm *r,sun_wferm *s)
*
*******************************************************************************/

#define DIRAC_WIL_C

#include"modules.h"



void applyWilsonDiracOperator(sun_wferm *r,sun_wferm *s)
{
#if(DIM==4)
   int n;
   int n0,nm0,n1,nm1,n2,nm2,n3,nm3;
   double fm;
   sun_wferm z1,z2,z3;

   fm=(4.+runParams.mass);

   for(n=0;n<VOL;n++)
   {
      sunwferm_real_mult(r[n],fm,s[n]);

      n0=neib[n][0];
      n1=neib[n][1];
      n2=neib[n][2];
      n3=neib[n][3];
      nm0=neib[n][DIM];
      nm1=neib[n][DIM+1];
      nm2=neib[n][DIM+2];
      nm3=neib[n][DIM+3];

      if(n0<n)
      {
         sunwferm_sun_mult(z1,*pu[n][0],s[n0]);
         mul_sunwferm_mg0(z2,z1);
         sunwferm_add_single(z1,z2);
         sunwferm_real_mult(z3,-1.,z1);
      }
      else
      {
         sunwferm_sun_mult(z3,*pu[n][0],s[n0]);
         mul_sunwferm_mg0(z2,z3);
         sunwferm_add_single(z3,z2);
      }

      if(nm0>n)
      {
         sunwferm_sun_dag_mult(z1,*pu[nm0][0],s[nm0]);
         mul_sunwferm_g0(z2,z1);
         sunwferm_add_single(z1,z2);
         sunwferm_real_mult(z2,-1.,z1);
         sunwferm_add_single(z3,z2);
      }
      else
      {
         sunwferm_sun_dag_mult(z1,*pu[nm0][0],s[nm0]);
         sunwferm_add_single(z3,z1);
         mul_sunwferm_g0(z2,z1);
         sunwferm_add_single(z3,z2);
      }

      sunwferm_sun_mult(z1,*pu[n][1],s[n1]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg1(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm1][1],s[nm1]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g1(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_mult(z1,*pu[n][2],s[n2]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg2(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm2][2],s[nm2]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g2(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_mult(z1,*pu[n][3],s[n3]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg3(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm3][3],s[nm3]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g3(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_real_mult(z1,-0.5,z3);
      sunwferm_add_single(r[n],z1);
   }
#else
   error(1,"apply_dirac_wil [dirac_wil.c]","DIM=4 is mandatory!");
#endif
}


void applyComplexConjDiracOperator(sun_wferm *r,sun_wferm *s)
{
#if(DIM==4)
   int n;
   int n0,nm0,n1,nm1,n2,nm2,n3,nm3;
   double fm;
   sun_wferm z1,z2,z3;

   fm=(4.+runParams.mass);

   for(n=0;n<VOL;n++)
   {
      sunwferm_real_mult(r[n],fm,s[n]);

      n0=neib[n][0];
      n1=neib[n][1];
      n2=neib[n][2];
      n3=neib[n][3];
      nm0=neib[n][DIM];
      nm1=neib[n][DIM+1];
      nm2=neib[n][DIM+2];
      nm3=neib[n][DIM+3];

      if(n0<n)
      {
         sunwferm_sun_mult(z1,*pu[n][0],s[n0]);
         mul_sunwferm_g0(z2,z1);
         sunwferm_add_single(z1,z2);
        sunwferm_real_mult(z3,-1.,z1);
      }
      else
      {
         sunwferm_sun_mult(z3,*pu[n][0],s[n0]);
         mul_sunwferm_g0(z2,z3);
         sunwferm_add_single(z3,z2);
      }

      if(nm0>n)
      {
         sunwferm_sun_dag_mult(z1,*pu[nm0][0],s[nm0]);
         mul_sunwferm_mg0(z2,z1);
         sunwferm_add_single(z1,z2);
        sunwferm_real_mult(z2,-1.,z1);
         sunwferm_add_single(z3,z2);
      }
      else
      {
         sunwferm_sun_dag_mult(z1,*pu[nm0][0],s[nm0]);
         sunwferm_add_single(z3,z1);
         mul_sunwferm_mg0(z2,z1);
         sunwferm_add_single(z3,z2);
      }

      sunwferm_sun_mult(z1,*pu[n][1],s[n1]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g1(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm1][1],s[nm1]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg1(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_mult(z1,*pu[n][2],s[n2]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g2(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm2][2],s[nm2]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg2(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_mult(z1,*pu[n][3],s[n3]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_g3(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_sun_dag_mult(z1,*pu[nm3][3],s[nm3]);
      sunwferm_add_single(z3,z1);
      mul_sunwferm_mg3(z2,z1);
      sunwferm_add_single(z3,z2);

      sunwferm_real_mult(z1,-0.5,z3);
      sunwferm_add_single(r[n],z1);
   }
#else
   error(1,"apply_dirac_wil [dirac_wil.c]","DIM=4 is mandatory!");
#endif
}
