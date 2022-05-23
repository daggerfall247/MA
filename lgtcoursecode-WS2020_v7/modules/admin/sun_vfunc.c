
/*******************************************************************************
*
* File sun_vfunc.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes routines to reunitarise the gauge field
* 
* The SU(3) routines included are similar to those used by Martin Luescher in
* the DD-HMC code.
* 
* Externally accessible functions:
* 
* void project_to_su3(su3mat *u)
*      Projects an approximate SU(3) matrix back to SU(3).
*
* void project_to_su2(su2mat *u)
*      Projects an approximate SU(2) matrix back to SU(2).
*
* void project_gfield_to_sun(sun_mat *u[VOL][DIM])
*      Projects approximate SU(2,3) matrices back to SU(2,3),
*      looping over all matrices in the gaugefield.
*
*******************************************************************************/

#define SUN_VFUNC_C

#include"modules.h"


#if(SUN==3)
static void normalize(su3vec *v)
{
   int i;
   double *r,fact;

   r=(double*)(v);
   fact=0.0;
   for(i=0;i<6;i++)
      fact+=r[i]*r[i];
   fact=1.0/sqrt(fact);
   for(i=0;i<6;i++)
      r[i]*=fact;
}
#endif


void project_to_su3(su3mat *u)
{
#if(SUN==3)
   su3vec *v1,*v2,*v3;

   v1=(su3vec*)(u);
   v2=v1+1;
   v3=v1+2;

   normalize(v1);
   su3vec_cross_prod(*v3,*v1,*v2);
   normalize(v3);
   su3vec_cross_prod(*v2,*v3,*v1);
#else
   error(1,"project_to_su3","Not in SU(3) mode!");
#endif
}


void project_to_su2(su2mat *u)
{
#if(SUN==2)
   double zw;
   
   su2_det_sqrt(zw,*u);
   su2_dble_div(*u,zw);
#else
   error(1,"project_to_su2","Not in SU(2) mode!");
#endif
}


void project_gfield_to_sun(sun_mat *u[VOL][DIM])
{
   int ii,jj;

   for(ii=0;ii<VOL;ii++)
      for(jj=0;jj<DIM;jj++)
      {
#if(SUN==2)
         project_to_su2(u[ii][jj]);
#elif(SUN==3)
         project_to_su3(u[ii][jj]);
#endif
      }
}
