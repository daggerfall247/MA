
/*******************************************************************************
*
* File gauge.h
*
* Copyright (C) 2013 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Includes definitions and macros for the gauge field variables
*
*******************************************************************************/

#define GAUGE_H

#ifndef COMPLEX_H
 #include"complex.h"
#endif

typedef struct
{
   double c0,c1,c2,c3;
} su2mat;

typedef struct
{
   double c1,c2,c3;
} su2alg;

typedef struct
{
   complex c11,c12,c13;
   complex c21,c22,c23;
   complex c31,c32,c33;
} su3mat;

typedef struct
{
   double c1,c2,c3,c4,c5,c6,c7,c8;
} su3alg;

typedef struct
{
   complex c1,c2,c3;
} su3vec;

/* SU(2) macros */

/* u=u^+ */
#define su2_dag_single(u){ \
        (u).c1=-(u).c1; \
        (u).c2=-(u).c2; \
        (u).c3=-(u).c3;}

/* u=v^+ */
#define su2_dag(u,v){ \
        (u).c0=(v).c0; \
        (u).c1=-(v).c1; \
        (u).c2=-(v).c2; \
        (u).c3=-(v).c3;}

/* u=u^* */
#define su2_star(u){ \
        (u).c1=-(u).c1; \
        (u).c3=-(u).c3;}

/* u=v^* */
#define su2_star_2g(u,v){ \
        (u).c0=(v).c0; \
        (u).c1=-(v).c1; \
        (u).c2=(v).c2; \
        (u).c3=-(v).c3;}

/* u=v^T */
#define su2_transp(u,v){ \
        (u).c0=(v).c0; \
        (u).c1=(v).c1; \
        (u).c2=-(v).c2; \
        (u).c3=(v).c3;}

/* u=v*w */
#define su2_mat_mul(u,v,w){ \
        (u).c0= (v).c0*(w).c0-(v).c1*(w).c1 \
               -(v).c2*(w).c2-(v).c3*(w).c3; \
        (u).c1= (v).c0*(w).c1+(v).c1*(w).c0 \
               -(v).c2*(w).c3+(v).c3*(w).c2; \
        (u).c2= (v).c0*(w).c2+(v).c1*(w).c3 \
               +(v).c2*(w).c0-(v).c3*(w).c1; \
        (u).c3= (v).c0*(w).c3-(v).c1*(w).c2 \
               +(v).c2*(w).c1+(v).c3*(w).c0;}

/* u=v*dag(w) */
#define su2_mat_mul_dag(u,v,w){ \
        (u).c0= (v).c0*(w).c0+(v).c1*(w).c1 \
               +(v).c2*(w).c2+(v).c3*(w).c3; \
        (u).c1=-(v).c0*(w).c1+(v).c1*(w).c0 \
               +(v).c2*(w).c3-(v).c3*(w).c2; \
        (u).c2=-(v).c0*(w).c2-(v).c1*(w).c3 \
               +(v).c2*(w).c0+(v).c3*(w).c1; \
        (u).c3=-(v).c0*(w).c3+(v).c1*(w).c2 \
               -(v).c2*(w).c1+(v).c3*(w).c0;}

/* u=0 */
#define su2_zero(u){ \
        (u).c0=0.; \
        (u).c1=0.; \
        (u).c2=0.; \
        (u).c3=0.;}

/* u=1 */
#define su2_unit(u){ \
        (u).c0=1.; \
        (u).c1=0.; \
        (u).c2=0.; \
        (u).c3=0.;}

/* u=v+w */
#define su2_add(u,v,w){ \
        (u).c0=(v).c0+(w).c0; \
        (u).c1=(v).c1+(w).c1; \
        (u).c2=(v).c2+(w).c2; \
        (u).c3=(v).c3+(w).c3;}

/* u=u+v */
#define su2_self_add(u,v){ \
        (u).c0+=(v).c0; \
        (u).c1+=(v).c1; \
        (u).c2+=(v).c2; \
        (u).c3+=(v).c3;}

/* u=v-w */
#define su2_sub(u,v,w){ \
        (u).c0=(v).c0-(w).c0; \
        (u).c1=(v).c1-(w).c1; \
        (u).c2=(v).c2-(w).c2; \
        (u).c3=(v).c3-(w).c3;}

/* u=u-v */
#define su2_self_sub(u,v){ \
        (u).c0-=(v).c0; \
        (u).c1-=(v).c1; \
        (u).c2-=(v).c2; \
        (u).c3-=(v).c3;}

/* u=sqrt(det(v)) */
#define su2_det_sqrt(u,v){ \
        (u)=(v).c0*(v).c0; \
        (u)+=(v).c1*(v).c1; \
        (u)+=(v).c2*(v).c2; \
        (u)+=(v).c3*(v).c3; \
        (u)=sqrt(fabs((u)));}

/* u=|det(v)| */
#define su2_det(u,v){ \
        (u)=(v).c0*(v).c0; \
        (u)+=(v).c1*(v).c1; \
        (u)+=(v).c2*(v).c2; \
        (u)+=(v).c3*(v).c3; \
        (u)=fabs(u);}

/* u=u*v : v real */
#define su2_dble_mul(u,v){ \
        (u).c0*=(v); \
        (u).c1*=(v); \
        (u).c2*=(v); \
        (u).c3*=(v);}

/* u=u/v : v real */
#define su2_dble_div(u,v){ \
        (u).c0/=(v); \
        (u).c1/=(v); \
        (u).c2/=(v); \
        (u).c3/=(v);}

/* tr=Re(Tr(u)) : tr real */
#define su2_trace(tr,u){ \
        (tr)=2.*(u).c0;}

/* tr=Tr(u) : tr complex */
#define su2_trace_cl(tr,u){ \
        (tr).re=2.*(u).c0; \
        (tr).im=0.;}

/* a=u : a complex 2x2 array */
#define mksu2_2x2(a,u){ \
        ((a)[0][0]).re=(u).c0; \
        ((a)[0][0]).im=(u).c3; \
        ((a)[0][1]).re=(u).c2; \
        ((a)[0][1]).im=(u).c1; \
        ((a)[1][0]).re=-(u).c2; \
        ((a)[1][0]).im=(u).c1; \
        ((a)[1][1]).re=(u).c0; \
        ((a)[1][1]).im=-(u).c3;}

/* a=u^+ : a complex 2x2 array */
#define mksu2_2x2_dag(a,u){ \
        ((a)[0][0]).re=(u).c0; \
        ((a)[0][0]).im=-(u).c3; \
        ((a)[0][1]).re=-(u).c2; \
        ((a)[0][1]).im=-(u).c1; \
        ((a)[1][0]).re=(u).c2; \
        ((a)[1][0]).im=-(u).c1; \
        ((a)[1][1]).re=(u).c0; \
        ((a)[1][1]).im=(u).c3;}

/* a=u : a double 4 array */
#define mk_su2_dble_array(a,u){ \
        *(a)=(u).c0; \
        *(a+1)=(u).c1; \
        *(a+2)=(u).c2; \
        *(a+3)=(u).c3;}

/* u=a : a double 4 array */
#define mk_dble_array_su2(a,u){ \
        (u).c0=*(a); \
        (u).c1=*(a+1); \
        (u).c2=*(a+2); \
        (u).c3=*(a+3);}
        
/* su(2) macros */

/* u=0 */
#define set_alg_zero_su2(u){ \
        (u).c1=0.; \
        (u).c2=0.; \
        (u).c3=0.;}

/* u=a : a double 3 array */
#define mk_dble_array_su2_alg(a,u){ \
        (u).c1=*(a); \
        (u).c2=*(a+1); \
        (u).c3=*(a+2);}

/* a=u : a double 3 array */
#define mk_su2_alg_dble_array(a,u){ \
        *(a)=(u).c1; \
        *(a+1)=(u).c2; \
        *(a+2)=(u).c3;}

/* r=r+c*s : c real */
#define su2_alg_mul_real_add_single(r,c,s){ \
        (r).c1+=(c)*(s).c1; \
        (r).c2+=(c)*(s).c2; \
        (r).c3+=(c)*(s).c3;}

/* SU(3) macros */

/* u=v^+ */
#define su3_dag(u,v){ \
        (u).c11.re= (v).c11.re; \
        (u).c11.im=-(v).c11.im; \
        (u).c12.re= (v).c21.re; \
        (u).c12.im=-(v).c21.im; \
        (u).c13.re= (v).c31.re; \
        (u).c13.im=-(v).c31.im; \
        (u).c21.re= (v).c12.re; \
        (u).c21.im=-(v).c12.im; \
        (u).c22.re= (v).c22.re; \
        (u).c22.im=-(v).c22.im; \
        (u).c23.re= (v).c32.re; \
        (u).c23.im=-(v).c32.im; \
        (u).c31.re= (v).c13.re; \
        (u).c31.im=-(v).c13.im; \
        (u).c32.re= (v).c23.re; \
        (u).c32.im=-(v).c23.im; \
        (u).c33.re= (v).c33.re; \
        (u).c33.im=-(v).c33.im;}

/* u=u^* */
#define su3_star(u){ \
        (u).c11.im=-(u).c11.im; \
        (u).c12.im=-(u).c12.im; \
        (u).c13.im=-(u).c13.im; \
        (u).c21.im=-(u).c21.im; \
        (u).c22.im=-(u).c22.im; \
        (u).c23.im=-(u).c23.im; \
        (u).c31.im=-(u).c31.im; \
        (u).c32.im=-(u).c32.im; \
        (u).c33.im=-(u).c33.im;}

/* u=v^* */
#define su3_star_2g(u,v){ \
        (u).c11.re= (v).c11.re; \
        (u).c11.im=-(v).c11.im; \
        (u).c12.re= (v).c12.re; \
        (u).c12.im=-(v).c12.im; \
        (u).c13.re= (v).c13.re; \
        (u).c13.im=-(v).c13.im; \
        (u).c21.re= (v).c21.re; \
        (u).c21.im=-(v).c21.im; \
        (u).c22.re= (v).c22.re; \
        (u).c22.im=-(v).c22.im; \
        (u).c23.re= (v).c23.re; \
        (u).c23.im=-(v).c23.im; \
        (u).c31.re= (v).c31.re; \
        (u).c31.im=-(v).c31.im; \
        (u).c32.re= (v).c32.re; \
        (u).c32.im=-(v).c32.im; \
        (u).c33.re= (v).c33.re; \
        (u).c33.im=-(v).c33.im;}

/* u=v^T */
#define su3_transp(u,v){ \
        (u).c11.re=(v).c11.re; \
        (u).c11.im=(v).c11.im; \
        (u).c12.re=(v).c21.re; \
        (u).c12.im=(v).c21.im; \
        (u).c13.re=(v).c31.re; \
        (u).c13.im=(v).c31.im; \
        (u).c21.re=(v).c12.re; \
        (u).c21.im=(v).c12.im; \
        (u).c22.re=(v).c22.re; \
        (u).c22.im=(v).c22.im; \
        (u).c23.re=(v).c32.re; \
        (u).c23.im=(v).c32.im; \
        (u).c31.re=(v).c13.re; \
        (u).c31.im=(v).c13.im; \
        (u).c32.re=(v).c23.re; \
        (u).c32.im=(v).c23.im; \
        (u).c33.re=(v).c33.re; \
        (u).c33.im=(v).c33.im;}

/* u=v*w */
#define su3_mat_mul(u,v,w){ \
        (u).c11.re= (v).c11.re*(w).c11.re-(v).c11.im*(w).c11.im  \
                   +(v).c12.re*(w).c21.re-(v).c12.im*(w).c21.im  \
                   +(v).c13.re*(w).c31.re-(v).c13.im*(w).c31.im; \
        (u).c11.im= (v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
                   +(v).c12.re*(w).c21.im+(v).c12.im*(w).c21.re  \
                   +(v).c13.re*(w).c31.im+(v).c13.im*(w).c31.re; \
        (u).c12.re= (v).c11.re*(w).c12.re-(v).c11.im*(w).c12.im  \
                   +(v).c12.re*(w).c22.re-(v).c12.im*(w).c22.im  \
                   +(v).c13.re*(w).c32.re-(v).c13.im*(w).c32.im; \
        (u).c12.im= (v).c11.re*(w).c12.im+(v).c11.im*(w).c12.re  \
                   +(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
                   +(v).c13.re*(w).c32.im+(v).c13.im*(w).c32.re; \
        (u).c13.re= (v).c11.re*(w).c13.re-(v).c11.im*(w).c13.im  \
                   +(v).c12.re*(w).c23.re-(v).c12.im*(w).c23.im  \
                   +(v).c13.re*(w).c33.re-(v).c13.im*(w).c33.im; \
        (u).c13.im= (v).c11.re*(w).c13.im+(v).c11.im*(w).c13.re  \
                   +(v).c12.re*(w).c23.im+(v).c12.im*(w).c23.re  \
                   +(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
        (u).c21.re= (v).c21.re*(w).c11.re-(v).c21.im*(w).c11.im  \
                   +(v).c22.re*(w).c21.re-(v).c22.im*(w).c21.im  \
                   +(v).c23.re*(w).c31.re-(v).c23.im*(w).c31.im; \
        (u).c21.im= (v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
                   +(v).c22.re*(w).c21.im+(v).c22.im*(w).c21.re  \
                   +(v).c23.re*(w).c31.im+(v).c23.im*(w).c31.re; \
        (u).c22.re= (v).c21.re*(w).c12.re-(v).c21.im*(w).c12.im  \
                   +(v).c22.re*(w).c22.re-(v).c22.im*(w).c22.im  \
                   +(v).c23.re*(w).c32.re-(v).c23.im*(w).c32.im; \
        (u).c22.im= (v).c21.re*(w).c12.im+(v).c21.im*(w).c12.re  \
                   +(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
                   +(v).c23.re*(w).c32.im+(v).c23.im*(w).c32.re; \
        (u).c23.re= (v).c21.re*(w).c13.re-(v).c21.im*(w).c13.im  \
                   +(v).c22.re*(w).c23.re-(v).c22.im*(w).c23.im  \
                   +(v).c23.re*(w).c33.re-(v).c23.im*(w).c33.im; \
        (u).c23.im= (v).c21.re*(w).c13.im+(v).c21.im*(w).c13.re  \
                   +(v).c22.re*(w).c23.im+(v).c22.im*(w).c23.re  \
                   +(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
        (u).c31.re= (v).c31.re*(w).c11.re-(v).c31.im*(w).c11.im  \
                   +(v).c32.re*(w).c21.re-(v).c32.im*(w).c21.im  \
                   +(v).c33.re*(w).c31.re-(v).c33.im*(w).c31.im; \
        (u).c31.im= (v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
                   +(v).c32.re*(w).c21.im+(v).c32.im*(w).c21.re  \
                   +(v).c33.re*(w).c31.im+(v).c33.im*(w).c31.re; \
        (u).c32.re= (v).c31.re*(w).c12.re-(v).c31.im*(w).c12.im  \
                   +(v).c32.re*(w).c22.re-(v).c32.im*(w).c22.im  \
                   +(v).c33.re*(w).c32.re-(v).c33.im*(w).c32.im; \
        (u).c32.im= (v).c31.re*(w).c12.im+(v).c31.im*(w).c12.re  \
                   +(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
                   +(v).c33.re*(w).c32.im+(v).c33.im*(w).c32.re; \
        (u).c33.re= (v).c31.re*(w).c13.re-(v).c31.im*(w).c13.im  \
                   +(v).c32.re*(w).c23.re-(v).c32.im*(w).c23.im  \
                   +(v).c33.re*(w).c33.re-(v).c33.im*(w).c33.im; \
        (u).c33.im= (v).c31.re*(w).c13.im+(v).c31.im*(w).c13.re  \
                   +(v).c32.re*(w).c23.im+(v).c32.im*(w).c23.re  \
                   +(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re;}

/* u=v*w^+ */
#define su3_mat_mul_dag(u,v,w){ \
        (u).c11.re= (v).c11.re*(w).c11.re+(v).c11.im*(w).c11.im  \
                   +(v).c12.re*(w).c12.re+(v).c12.im*(w).c12.im  \
                   +(v).c13.re*(w).c13.re+(v).c13.im*(w).c13.im; \
        (u).c11.im=-(v).c11.re*(w).c11.im+(v).c11.im*(w).c11.re  \
                   -(v).c12.re*(w).c12.im+(v).c12.im*(w).c12.re  \
                   -(v).c13.re*(w).c13.im+(v).c13.im*(w).c13.re; \
        (u).c12.re= (v).c11.re*(w).c21.re+(v).c11.im*(w).c21.im  \
                   +(v).c12.re*(w).c22.re+(v).c12.im*(w).c22.im  \
                   +(v).c13.re*(w).c23.re+(v).c13.im*(w).c23.im; \
        (u).c12.im=-(v).c11.re*(w).c21.im+(v).c11.im*(w).c21.re  \
                   -(v).c12.re*(w).c22.im+(v).c12.im*(w).c22.re  \
                   -(v).c13.re*(w).c23.im+(v).c13.im*(w).c23.re; \
        (u).c13.re= (v).c11.re*(w).c31.re+(v).c11.im*(w).c31.im  \
                   +(v).c12.re*(w).c32.re+(v).c12.im*(w).c32.im  \
                   +(v).c13.re*(w).c33.re+(v).c13.im*(w).c33.im; \
        (u).c13.im=-(v).c11.re*(w).c31.im+(v).c11.im*(w).c31.re  \
                   -(v).c12.re*(w).c32.im+(v).c12.im*(w).c32.re  \
                   -(v).c13.re*(w).c33.im+(v).c13.im*(w).c33.re; \
        (u).c21.re= (v).c21.re*(w).c11.re+(v).c21.im*(w).c11.im  \
                   +(v).c22.re*(w).c12.re+(v).c22.im*(w).c12.im  \
                   +(v).c23.re*(w).c13.re+(v).c23.im*(w).c13.im; \
        (u).c21.im=-(v).c21.re*(w).c11.im+(v).c21.im*(w).c11.re  \
                   -(v).c22.re*(w).c12.im+(v).c22.im*(w).c12.re  \
                   -(v).c23.re*(w).c13.im+(v).c23.im*(w).c13.re; \
        (u).c22.re= (v).c21.re*(w).c21.re+(v).c21.im*(w).c21.im  \
                   +(v).c22.re*(w).c22.re+(v).c22.im*(w).c22.im  \
                   +(v).c23.re*(w).c23.re+(v).c23.im*(w).c23.im; \
        (u).c22.im=-(v).c21.re*(w).c21.im+(v).c21.im*(w).c21.re  \
                   -(v).c22.re*(w).c22.im+(v).c22.im*(w).c22.re  \
                   -(v).c23.re*(w).c23.im+(v).c23.im*(w).c23.re; \
        (u).c23.re= (v).c21.re*(w).c31.re+(v).c21.im*(w).c31.im  \
                   +(v).c22.re*(w).c32.re+(v).c22.im*(w).c32.im  \
                   +(v).c23.re*(w).c33.re+(v).c23.im*(w).c33.im; \
        (u).c23.im=-(v).c21.re*(w).c31.im+(v).c21.im*(w).c31.re  \
                   -(v).c22.re*(w).c32.im+(v).c22.im*(w).c32.re  \
                   -(v).c23.re*(w).c33.im+(v).c23.im*(w).c33.re; \
        (u).c31.re= (v).c31.re*(w).c11.re+(v).c31.im*(w).c11.im  \
                   +(v).c32.re*(w).c12.re+(v).c32.im*(w).c12.im  \
                   +(v).c33.re*(w).c13.re+(v).c33.im*(w).c13.im; \
        (u).c31.im=-(v).c31.re*(w).c11.im+(v).c31.im*(w).c11.re  \
                   -(v).c32.re*(w).c12.im+(v).c32.im*(w).c12.re  \
                   -(v).c33.re*(w).c13.im+(v).c33.im*(w).c13.re; \
        (u).c32.re= (v).c31.re*(w).c21.re+(v).c31.im*(w).c21.im  \
                   +(v).c32.re*(w).c22.re+(v).c32.im*(w).c22.im  \
                   +(v).c33.re*(w).c23.re+(v).c33.im*(w).c23.im; \
        (u).c32.im=-(v).c31.re*(w).c21.im+(v).c31.im*(w).c21.re  \
                   -(v).c32.re*(w).c22.im+(v).c32.im*(w).c22.re  \
                   -(v).c33.re*(w).c23.im+(v).c33.im*(w).c23.re; \
        (u).c33.re= (v).c31.re*(w).c31.re+(v).c31.im*(w).c31.im  \
                   +(v).c32.re*(w).c32.re+(v).c32.im*(w).c32.im  \
                   +(v).c33.re*(w).c33.re+(v).c33.im*(w).c33.im; \
        (u).c33.im=-(v).c31.re*(w).c31.im+(v).c31.im*(w).c31.re  \
                   -(v).c32.re*(w).c32.im+(v).c32.im*(w).c32.re  \
                   -(v).c33.re*(w).c33.im+(v).c33.im*(w).c33.re;}

/* u=0 */
#define su3_zero(u){ \
        (u).c11.re=0.; \
        (u).c11.im=0.; \
        (u).c12.re=0.; \
        (u).c12.im=0.; \
        (u).c13.re=0.; \
        (u).c13.im=0.; \
        (u).c21.re=0.; \
        (u).c21.im=0.; \
        (u).c22.re=0.; \
        (u).c22.im=0.; \
        (u).c23.re=0.; \
        (u).c23.im=0.; \
        (u).c31.re=0.; \
        (u).c31.im=0.; \
        (u).c32.re=0.; \
        (u).c32.im=0.; \
        (u).c33.re=0.; \
        (u).c33.im=0.;}

/* u=1 */
#define su3_unit(u){ \
        (u).c11.re=1.; \
        (u).c11.im=0.; \
        (u).c12.re=0.; \
        (u).c12.im=0.; \
        (u).c13.re=0.; \
        (u).c13.im=0.; \
        (u).c21.re=0.; \
        (u).c21.im=0.; \
        (u).c22.re=1.; \
        (u).c22.im=0.; \
        (u).c23.re=0.; \
        (u).c23.im=0.; \
        (u).c31.re=0.; \
        (u).c31.im=0.; \
        (u).c32.re=0.; \
        (u).c32.im=0.; \
        (u).c33.re=1.; \
        (u).c33.im=0.;}

/* u=v+w */
#define su3_add(u,v,w){ \
        (u).c11.re=(v).c11.re+(w).c11.re; \
        (u).c11.im=(v).c11.im+(w).c11.im; \
        (u).c12.re=(v).c12.re+(w).c12.re; \
        (u).c12.im=(v).c12.im+(w).c12.im; \
        (u).c13.re=(v).c13.re+(w).c13.re; \
        (u).c13.im=(v).c13.im+(w).c13.im; \
        (u).c21.re=(v).c21.re+(w).c21.re; \
        (u).c21.im=(v).c21.im+(w).c21.im; \
        (u).c22.re=(v).c22.re+(w).c22.re; \
        (u).c22.im=(v).c22.im+(w).c22.im; \
        (u).c23.re=(v).c23.re+(w).c23.re; \
        (u).c23.im=(v).c23.im+(w).c23.im; \
        (u).c31.re=(v).c31.re+(w).c31.re; \
        (u).c31.im=(v).c31.im+(w).c31.im; \
        (u).c32.re=(v).c32.re+(w).c32.re; \
        (u).c32.im=(v).c32.im+(w).c32.im; \
        (u).c33.re=(v).c33.re+(w).c33.re; \
        (u).c33.im=(v).c33.im+(w).c33.im;}

/* u=u+v */
#define su3_self_add(u,v){ \
        (u).c11.re+=(v).c11.re; \
        (u).c11.im+=(v).c11.im; \
        (u).c12.re+=(v).c12.re; \
        (u).c12.im+=(v).c12.im; \
        (u).c13.re+=(v).c13.re; \
        (u).c13.im+=(v).c13.im; \
        (u).c21.re+=(v).c21.re; \
        (u).c21.im+=(v).c21.im; \
        (u).c22.re+=(v).c22.re; \
        (u).c22.im+=(v).c22.im; \
        (u).c23.re+=(v).c23.re; \
        (u).c23.im+=(v).c23.im; \
        (u).c31.re+=(v).c31.re; \
        (u).c31.im+=(v).c31.im; \
        (u).c32.re+=(v).c32.re; \
        (u).c32.im+=(v).c32.im; \
        (u).c33.re+=(v).c33.re; \
        (u).c33.im+=(v).c33.im;}

/* u=v-w */
#define su3_sub(u,v,w){ \
        (u).c11.re=(v).c11.re-(w).c11.re; \
        (u).c11.im=(v).c11.im-(w).c11.im; \
        (u).c12.re=(v).c12.re-(w).c12.re; \
        (u).c12.im=(v).c12.im-(w).c12.im; \
        (u).c13.re=(v).c13.re-(w).c13.re; \
        (u).c13.im=(v).c13.im-(w).c13.im; \
        (u).c21.re=(v).c21.re-(w).c21.re; \
        (u).c21.im=(v).c21.im-(w).c21.im; \
        (u).c22.re=(v).c22.re-(w).c22.re; \
        (u).c22.im=(v).c22.im-(w).c22.im; \
        (u).c23.re=(v).c23.re-(w).c23.re; \
        (u).c23.im=(v).c23.im-(w).c23.im; \
        (u).c31.re=(v).c31.re-(w).c31.re; \
        (u).c31.im=(v).c31.im-(w).c31.im; \
        (u).c32.re=(v).c32.re-(w).c32.re; \
        (u).c32.im=(v).c32.im-(w).c32.im; \
        (u).c33.re=(v).c33.re-(w).c33.re; \
        (u).c33.im=(v).c33.im-(w).c33.im;}

/* u=u-v */
#define su3_self_sub(u,v){ \
        (u).c11.re+=-(v).c11.re; \
        (u).c11.im+=-(v).c11.im; \
        (u).c12.re+=-(v).c12.re; \
        (u).c12.im+=-(v).c12.im; \
        (u).c13.re+=-(v).c13.re; \
        (u).c13.im+=-(v).c13.im; \
        (u).c21.re+=-(v).c21.re; \
        (u).c21.im+=-(v).c21.im; \
        (u).c22.re+=-(v).c22.re; \
        (u).c22.im+=-(v).c22.im; \
        (u).c23.re+=-(v).c23.re; \
        (u).c23.im+=-(v).c23.im; \
        (u).c31.re+=-(v).c31.re; \
        (u).c31.im+=-(v).c31.im; \
        (u).c32.re+=-(v).c32.re; \
        (u).c32.im+=-(v).c32.im; \
        (u).c33.re+=-(v).c33.re; \
        (u).c33.im+=-(v).c33.im;}

/* u=u*z : z real */
#define su3_dble_mul(u,z){ \
        (u).c11.re*=(z); \
        (u).c11.im*=(z); \
        (u).c12.re*=(z); \
        (u).c12.im*=(z); \
        (u).c13.re*=(z); \
        (u).c13.im*=(z); \
        (u).c21.re*=(z); \
        (u).c21.im*=(z); \
        (u).c22.re*=(z); \
        (u).c22.im*=(z); \
        (u).c23.re*=(z); \
        (u).c23.im*=(z); \
        (u).c31.re*=(z); \
        (u).c31.im*=(z); \
        (u).c32.re*=(z); \
        (u).c32.im*=(z); \
        (u).c33.re*=(z); \
        (u).c33.im*=(z);}

/* u=u/z : z real */
#define su3_dble_div(u,z){ \
        (u).c11.re/=(z); \
        (u).c11.im/=(z); \
        (u).c12.re/=(z); \
        (u).c12.im/=(z); \
        (u).c13.re/=(z); \
        (u).c13.im/=(z); \
        (u).c21.re/=(z); \
        (u).c21.im/=(z); \
        (u).c22.re/=(z); \
        (u).c22.im/=(z); \
        (u).c23.re/=(z); \
        (u).c23.im/=(z); \
        (u).c31.re/=(z); \
        (u).c31.im/=(z); \
        (u).c32.re/=(z); \
        (u).c32.im/=(z); \
        (u).c33.re/=(z); \
        (u).c33.im/=(z);}

/* u=v-w */
#define su3_sub(u,v,w){ \
        (u).c11.re=(v).c11.re-(w).c11.re; \
        (u).c11.im=(v).c11.im-(w).c11.im; \
        (u).c12.re=(v).c12.re-(w).c12.re; \
        (u).c12.im=(v).c12.im-(w).c12.im; \
        (u).c13.re=(v).c13.re-(w).c13.re; \
        (u).c13.im=(v).c13.im-(w).c13.im; \
        (u).c21.re=(v).c21.re-(w).c21.re; \
        (u).c21.im=(v).c21.im-(w).c21.im; \
        (u).c22.re=(v).c22.re-(w).c22.re; \
        (u).c22.im=(v).c22.im-(w).c22.im; \
        (u).c23.re=(v).c23.re-(w).c23.re; \
        (u).c23.im=(v).c23.im-(w).c23.im; \
        (u).c31.re=(v).c31.re-(w).c31.re; \
        (u).c31.im=(v).c31.im-(w).c31.im; \
        (u).c32.re=(v).c32.re-(w).c32.re; \
        (u).c32.im=(v).c32.im-(w).c32.im; \
        (u).c33.re=(v).c33.re-(w).c33.re; \
        (u).c33.im=(v).c33.im-(w).c33.im;}

/* tr=Re(Tr(u)) : tr real */
#define su3_trace_re(tr,u){ \
        (tr)=(u).c11.re+(u).c22.re+(u).c33.re;}

/* tr=Tr(u) : tr complex */
#define su3_trace(tr,u){ \
        (tr).re=(u).c11.re+(u).c22.re+(u).c33.re; \
        (tr).im=(u).c11.im+(u).c22.im+(u).c33.im;}

/* a=u : a complex 3x3 matrix */
#define mksu3_3x3(a,u){ \
        (a)[0][0]=(u).c11; \
        (a)[0][1]=(u).c12; \
        (a)[0][2]=(u).c13; \
        (a)[1][0]=(u).c21; \
        (a)[1][1]=(u).c22; \
        (a)[1][2]=(u).c23; \
        (a)[2][0]=(u).c31; \
        (a)[2][1]=(u).c32; \
        (a)[2][2]=(u).c33;}

/* a=u^+ : a complex 3x3 matrix */
#define mksu3_3x3_dag(a,u){ \
        (a)[0][0].re=(u).c11.re; \
        (a)[0][0].im=-(u).c11.im; \
        (a)[0][1].re=(u).c21.re; \
        (a)[0][1].im=-(u).c21.im; \
        (a)[0][2].re=(u).c31.re; \
        (a)[0][2].im=-(u).c31.im; \
        (a)[1][0].re=(u).c12.re; \
        (a)[1][0].im=-(u).c12.im; \
        (a)[1][1].re=(u).c22.re; \
        (a)[1][1].im=-(u).c22.im; \
        (a)[1][2].re=(u).c32.re; \
        (a)[1][2].im=-(u).c32.im; \
        (a)[2][0].re=(u).c13.re; \
        (a)[2][0].im=-(u).c13.im; \
        (a)[2][1].re=(u).c23.re; \
        (a)[2][1].im=-(u).c23.im; \
        (a)[2][2].re=(u).c33.re; \
        (a)[2][2].im=-(u).c33.im;}

/* a=u : a double 18 array */
#define mk_su3_dble_array(a,u){ \
        *(a)=(u).c11.re; \
        *(a+1)=(u).c11.im; \
        *(a+2)=(u).c12.re; \
        *(a+3)=(u).c12.im; \
        *(a+4)=(u).c13.re; \
        *(a+5)=(u).c13.im; \
        *(a+6)=(u).c21.re; \
        *(a+7)=(u).c21.im; \
        *(a+8)=(u).c22.re; \
        *(a+9)=(u).c22.im; \
        *(a+10)=(u).c23.re; \
        *(a+11)=(u).c23.im; \
        *(a+12)=(u).c31.re; \
        *(a+13)=(u).c31.im; \
        *(a+14)=(u).c32.re; \
        *(a+15)=(u).c32.im; \
        *(a+16)=(u).c33.re; \
        *(a+17)=(u).c33.im;}

/* u=a : a double 18 array */
#define mk_dble_array_su3(a,u){ \
        (u).c11.re=*(a); \
        (u).c11.im=*(a+1); \
        (u).c12.re=*(a+2); \
        (u).c12.im=*(a+3); \
        (u).c13.re=*(a+4); \
        (u).c13.im=*(a+5); \
        (u).c21.re=*(a+6); \
        (u).c21.im=*(a+7); \
        (u).c22.re=*(a+8); \
        (u).c22.im=*(a+9); \
        (u).c23.re=*(a+10); \
        (u).c23.im=*(a+11); \
        (u).c31.re=*(a+12); \
        (u).c31.im=*(a+13); \
        (u).c32.re=*(a+14); \
        (u).c32.im=*(a+15); \
        (u).c33.re=*(a+16); \
        (u).c33.im=*(a+17);}

/* su(3) macros */

/* u=0 */
#define set_alg_zero_su3(u){ \
        (u).c1=0.; \
        (u).c2=0.; \
        (u).c3=0.; \
        (u).c4=0.; \
        (u).c5=0.; \
        (u).c6=0.; \
        (u).c7=0.; \
        (u).c8=0.;}

/* u=a : a real 8 array */
#define mk_dble_array_su3_alg(a,u){ \
        (u).c1=*(a); \
        (u).c2=*(a+1); \
        (u).c3=*(a+2); \
        (u).c4=*(a+3); \
        (u).c5=*(a+4); \
        (u).c6=*(a+5); \
        (u).c7=*(a+6); \
        (u).c8=*(a+7);}

/* a=u : a real 8 array */
#define mk_su3_alg_dble_array(a,u){ \
        *(a)=(u).c1; \
        *(a+1)=(u).c2; \
        *(a+2)=(u).c3; \
        *(a+3)=(u).c4; \
        *(a+4)=(u).c5; \
        *(a+5)=(u).c6; \
        *(a+6)=(u).c7; \
        *(a+7)=(u).c8;}

/* r=r+a*s : a real */
#define su3_alg_mul_real_add_single(r,a,s){ \
        (r).c1+=(a)*(s).c1; \
        (r).c2+=(a)*(s).c2; \
        (r).c3+=(a)*(s).c3; \
        (r).c4+=(a)*(s).c4; \
        (r).c5+=(a)*(s).c5; \
        (r).c6+=(a)*(s).c6; \
        (r).c7+=(a)*(s).c7; \
        (r).c8+=(a)*(s).c8;}
        
/* Cross product of two SU(3) vectors fields -> https://math.stackexchange.com/a/129244 */
#define su3vec_cross_prod(v,w,z){ \
   (v).c1.re= (w).c2.re*(z).c3.re-(w).c2.im*(z).c3.im  \
             -(w).c3.re*(z).c2.re+(w).c3.im*(z).c2.im; \
   (v).c1.im= (w).c3.re*(z).c2.im+(w).c3.im*(z).c2.re  \
             -(w).c2.re*(z).c3.im-(w).c2.im*(z).c3.re; \
   (v).c2.re= (w).c3.re*(z).c1.re-(w).c3.im*(z).c1.im  \
             -(w).c1.re*(z).c3.re+(w).c1.im*(z).c3.im; \
   (v).c2.im= (w).c1.re*(z).c3.im+(w).c1.im*(z).c3.re  \
             -(w).c3.re*(z).c1.im-(w).c3.im*(z).c1.re; \
   (v).c3.re= (w).c1.re*(z).c2.re-(w).c1.im*(z).c2.im  \
             -(w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (v).c3.im= (w).c2.re*(z).c1.im+(w).c2.im*(z).c1.re  \
             -(w).c1.re*(z).c2.im-(w).c1.im*(z).c2.re;}

/*  Cabmar subgroups macros for SU(3) */
#define extr_sub1(u,v){ \
        (u).c0=((v).c11.re+(v).c22.re)/2.; \
        (u).c1=((v).c12.im+(v).c21.im)/2.; \
        (u).c2=((v).c12.re-(v).c21.re)/2.; \
        (u).c3=((v).c11.im-(v).c22.im)/2.;}

#define extr_sub2(u,v){ \
        (u).c0=((v).c22.re+(v).c33.re)/2.; \
        (u).c1=((v).c23.im+(v).c32.im)/2.; \
        (u).c2=((v).c23.re-(v).c32.re)/2.; \
        (u).c3=((v).c22.im-(v).c33.im)/2.;}

#define extr_sub3(u,v){ \
        (u).c0=((v).c11.re+(v).c33.re)/2.; \
        (u).c1=((v).c13.im+(v).c31.im)/2.; \
        (u).c2=((v).c13.re-(v).c31.re)/2.; \
        (u).c3=((v).c11.im-(v).c33.im)/2.;}

#define extr_sub1_dag(u,v){ \
        (u).c0=((v).c11.re+(v).c22.re)/2.; \
        (u).c1=-((v).c12.im+(v).c21.im)/2.; \
        (u).c2=-((v).c12.re-(v).c21.re)/2.; \
        (u).c3=-((v).c11.im-(v).c22.im)/2.;}

#define extr_sub2_dag(u,v){ \
        (u).c0=((v).c22.re+(v).c33.re)/2.; \
        (u).c1=-((v).c23.im+(v).c32.im)/2.; \
        (u).c2=-((v).c23.re-(v).c32.re)/2.; \
        (u).c3=-((v).c22.im-(v).c33.im)/2.;}

#define extr_sub3_dag(u,v){ \
        (u).c0=((v).c11.re+(v).c33.re)/2.; \
        (u).c1=-((v).c13.im+(v).c31.im)/2.; \
        (u).c2=-((v).c13.re-(v).c31.re)/2.; \
        (u).c3=-((v).c11.im-(v).c33.im)/2.;}

#define create_su3_sub1(u,v){ \
        (u).c11.re=(v).c0; \
        (u).c11.im=(v).c3; \
        (u).c12.re=(v).c2; \
        (u).c12.im=(v).c1; \
        (u).c13.re=0.; \
        (u).c13.im=0.; \
        (u).c21.re=-(v).c2; \
        (u).c21.im=(v).c1; \
        (u).c22.re=(v).c0; \
        (u).c22.im=-(v).c3; \
        (u).c23.re=0.; \
        (u).c23.im=0.; \
        (u).c31.re=0.; \
        (u).c31.im=0.; \
        (u).c32.re=0.; \
        (u).c32.im=0.; \
        (u).c33.re=1.; \
        (u).c33.im=0.;}

#define create_su3_sub2(u,v){ \
        (u).c11.re=1.; \
        (u).c11.im=0.; \
        (u).c12.re=0.; \
        (u).c12.im=0.; \
        (u).c13.re=0.; \
        (u).c13.im=0.; \
        (u).c21.re=0.; \
        (u).c21.im=0.; \
        (u).c22.re=(v).c0; \
        (u).c22.im=(v).c3; \
        (u).c23.re=(v).c2; \
        (u).c23.im=(v).c1; \
        (u).c31.re=0.; \
        (u).c31.im=0.; \
        (u).c32.re=-(v).c2; \
        (u).c32.im=(v).c1; \
        (u).c33.re=(v).c0; \
        (u).c33.im=-(v).c3;}

#define create_su3_sub3(u,v){ \
        (u).c11.re=(v).c0; \
        (u).c11.im=(v).c3; \
        (u).c12.re=0.; \
        (u).c12.im=0.; \
        (u).c13.re=(v).c2; \
        (u).c13.im=(v).c1; \
        (u).c21.re=0.; \
        (u).c21.im=0.; \
        (u).c22.re=1.; \
        (u).c22.im=0.; \
        (u).c23.re=0.; \
        (u).c23.im=0.; \
        (u).c31.re=-(v).c2; \
        (u).c31.im=(v).c1; \
        (u).c32.re=0.; \
        (u).c32.im=0.; \
        (u).c33.re=(v).c0; \
        (u).c33.im=-(v).c3;}
