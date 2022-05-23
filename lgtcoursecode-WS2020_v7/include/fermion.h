
/*******************************************************************************
*
* File fermion.h
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Includes definitions and macros for the pseudo fermion fields
*
*******************************************************************************/

#define FERMION_H

#ifndef COMPLEX_H
 #include"complex.h"
#endif

#ifndef GAUGE_H
 #include"gauge.h"
#endif

typedef struct
{
   complex c1,c2;
} su2vec;

typedef struct
{
   su2vec d1,d2,d3,d4;
} su2wferm;

typedef struct
{
   su3vec d1,d2,d3,d4;
} su3wferm;

/* SU(2) vector macros */

/* v=i*w */
#define su2_vec_mul_i(v,w){ \
        (v).c1.re=-(w).c1.im; \
        (v).c1.im=(w).c1.re; \
        (v).c2.re=-(w).c2.im; \
        (v).c2.im=(w).c2.re;}

/* v=(-i)*w */
#define su2_vec_mul_mi(v,w){ \
        (v).c1.re=(w).c1.im; \
        (v).c1.im=-(w).c1.re; \
        (v).c2.re=(w).c2.im; \
        (v).c2.im=-(w).c2.re;}

/* v=(-1)*w */
#define su2_vec_mul_m1(v,w){ \
        (v).c1.re=-(w).c1.re; \
        (v).c1.im=-(w).c1.im; \
        (v).c2.re=-(w).c2.re; \
        (v).c2.im=-(w).c2.im;}

/* r=s1+s2 */
#define su2vec_add(r,s1,s2){ \
        (r).c1.re=(s1).c1.re+(s2).c1.re; \
        (r).c1.im=(s1).c1.im+(s2).c1.im; \
        (r).c2.re=(s1).c2.re+(s2).c2.re; \
        (r).c2.im=(s1).c2.im+(s2).c2.im;}

/* r+=s */
#define su2vec_add_single(r,s){ \
        (r).c1.re+=(s).c1.re; \
        (r).c1.im+=(s).c1.im; \
        (r).c2.re+=(s).c2.re; \
        (r).c2.im+=(s).c2.im;}

/* r=s1-s2 */
#define su2vec_sub(r,s1,s2){ \
        (r).c1.re=(s1).c1.re-(s2).c1.re; \
        (r).c1.im=(s1).c1.im-(s2).c1.im; \
        (r).c2.re=(s1).c2.re-(s2).c2.re; \
        (r).c2.im=(s1).c2.im-(s2).c2.im;}

/* r-=s */
#define su2vec_sub_single(r,s){ \
        (r).c1.re-=(s).c1.re; \
        (r).c1.im-=(s).c1.im; \
        (r).c2.re-=(s).c2.re; \
        (r).c2.im-=(s).c2.im;}

/* r=c*s (c real) */
#define su2vec_real_mult(r,c,s){ \
        (r).c1.re=(c)*(s).c1.re; \
        (r).c1.im=(c)*(s).c1.im; \
        (r).c2.re=(c)*(s).c2.re; \
        (r).c2.im=(c)*(s).c2.im;}

/* r=u*s (u: SU(2) matrix; s: SU(2) vector) */
#define su2vec_su2_mult(r,u,s){ \
        (r).c1.re= (u).c0*(s).c1.re-(u).c3*(s).c1.im+(u).c2*(s).c2.re-(u).c1*(s).c2.im; \
        (r).c1.im= (u).c3*(s).c1.re+(u).c0*(s).c1.im+(u).c1*(s).c2.re+(u).c2*(s).c2.im; \
        (r).c2.re=-(u).c2*(s).c1.re-(u).c1*(s).c1.im+(u).c0*(s).c2.re+(u).c3*(s).c2.im; \
        (r).c2.im= (u).c1*(s).c1.re-(u).c2*(s).c1.im-(u).c3*(s).c2.re+(u).c0*(s).c2.im;}

/* r=u^+*s (u: SU(2) matrix; s: SU(2) vector) */
#define su2vec_su2_dag_mult(r,u,s){ \
        (r).c1.re= (u).c0*(s).c1.re+(u).c3*(s).c1.im-(u).c2*(s).c2.re+(u).c1*(s).c2.im; \
        (r).c1.im=-(u).c3*(s).c1.re+(u).c0*(s).c1.im-(u).c1*(s).c2.re-(u).c2*(s).c2.im; \
        (r).c2.re= (u).c2*(s).c1.re+(u).c1*(s).c1.im+(u).c0*(s).c2.re-(u).c3*(s).c2.im; \
        (r).c2.im=-(u).c1*(s).c1.re+(u).c2*(s).c1.im+(u).c3*(s).c2.re+(u).c0*(s).c2.im;}


/* SU(2) fermion macros */

/* r=s1+s2 */
#define su2wferm_add(r,s1,s2){ \
        su2vec_add((r).d1,(s1).d1,(s2).d1); \
        su2vec_add((r).d2,(s1).d2,(s2).d2); \
        su2vec_add((r).d3,(s1).d3,(s2).d3); \
        su2vec_add((r).d4,(s1).d4,(s2).d4);}

/* r+=s */
#define su2wferm_add_single(r,s){ \
        su2vec_add_single((r).d1,(s).d1); \
        su2vec_add_single((r).d2,(s).d2); \
        su2vec_add_single((r).d3,(s).d3); \
        su2vec_add_single((r).d4,(s).d4);}

/* r=s1-s2 */
#define su2wferm_sub(r,s1,s2){ \
        su2vec_sub((r).d1,(s1).d1,(s2).d1); \
        su2vec_sub((r).d2,(s1).d2,(s2).d2); \
        su2vec_sub((r).d3,(s1).d3,(s2).d3); \
        su2vec_sub((r).d4,(s1).d4,(s2).d4);}

/* r-=s */
#define su2wferm_sub_single(r,s){ \
        su2vec_sub_single((r).d1,(s).d1); \
        su2vec_sub_single((r).d2,(s).d2); \
        su2vec_sub_single((r).d3,(s).d3); \
        su2vec_sub_single((r).d4,(s).d4);}

/* r=c*s (c real) */
#define su2wferm_real_mult(r,c,s){ \
        su2vec_real_mult((r).d1,c,(s).d1); \
        su2vec_real_mult((r).d2,c,(s).d2); \
        su2vec_real_mult((r).d3,c,(s).d3); \
        su2vec_real_mult((r).d4,c,(s).d4);}

/* v=gamma_0*w */
#define mul_su2wferm_g0(v,w){ \
        su2_vec_mul_m1((v).d1,(w).d3); \
        su2_vec_mul_m1((v).d2,(w).d4); \
        su2_vec_mul_m1((v).d3,(w).d1); \
        su2_vec_mul_m1((v).d4,(w).d2);}

/* v=gamma_1*w */
#define mul_su2wferm_g1(v,w){ \
        su2_vec_mul_mi((v).d1,(w).d4); \
        su2_vec_mul_mi((v).d2,(w).d3); \
        su2_vec_mul_i((v).d3,(w).d2); \
        su2_vec_mul_i((v).d4,(w).d1);}

/* v=gamma_2*w */
#define mul_su2wferm_g2(v,w){ \
        su2_vec_mul_m1((v).d1,(w).d4); \
        (v).d2=(w).d3; \
        (v).d3=(w).d2; \
        su2_vec_mul_m1((v).d4,(w).d1);}

/* v=gamma_3*w */
#define mul_su2wferm_g3(v,w){ \
        su2_vec_mul_mi((v).d1,(w).d3); \
        su2_vec_mul_i((v).d2,(w).d4); \
        su2_vec_mul_i((v).d3,(w).d1); \
        su2_vec_mul_mi((v).d4,(w).d2);}

/* v=-gamma_0*w */
#define mul_su2wferm_mg0(v,w){ \
        (v).d1=(w).d3; \
        (v).d2=(w).d4; \
        (v).d3=(w).d1; \
        (v).d4=(w).d2;}

/* v=-gamma_1*w */
#define mul_su2wferm_mg1(v,w){ \
        su2_vec_mul_i((v).d1,(w).d4); \
        su2_vec_mul_i((v).d2,(w).d3); \
        su2_vec_mul_mi((v).d3,(w).d2); \
        su2_vec_mul_mi((v).d4,(w).d1);}

/* v=-gamma_2*w */
#define mul_su2wferm_mg2(v,w){ \
        (v).d1=(w).d4; \
        su2_vec_mul_m1((v).d2,(w).d3); \
        su2_vec_mul_m1((v).d3,(w).d2); \
        (v).d4=(w).d1;}

/* v=-gamma_3*w */
#define mul_su2wferm_mg3(v,w){ \
        su2_vec_mul_i((v).d1,(w).d3); \
        su2_vec_mul_mi((v).d2,(w).d4); \
        su2_vec_mul_mi((v).d3,(w).d1); \
        su2_vec_mul_i((v).d4,(w).d2);}

/* v=gamma_5*w */
#define mul_su2wferm_g5(v,w){ \
        (v).d1=(w).d1; \
        (v).d2=(w).d2; \
        su2_vec_mul_m1((v).d3,(w).d3); \
        su2_vec_mul_m1((v).d4,(w).d4);}

/* v=u*w (u: SU(2) matrix; w: fermion) */
#define su2wferm_su2_mult(v,u,w){ \
        su2vec_su2_mult((v).d1,u,(w).d1); \
        su2vec_su2_mult((v).d2,u,(w).d2); \
        su2vec_su2_mult((v).d3,u,(w).d3); \
        su2vec_su2_mult((v).d4,u,(w).d4);}

/* v=u^+*w (u: SU(2) matrix; w: fermion) */
#define su2wferm_su2_dag_mult(v,u,w){ \
        su2vec_su2_dag_mult((v).d1,u,(w).d1); \
        su2vec_su2_dag_mult((v).d2,u,(w).d2); \
        su2vec_su2_dag_mult((v).d3,u,(w).d3); \
        su2vec_su2_dag_mult((v).d4,u,(w).d4);}



/* SU(3) vector macros */

/* v=i*w */
#define su3_vec_mul_i(v,w){ \
        (v).c1.re=-(w).c1.im; \
        (v).c1.im=(w).c1.re; \
        (v).c2.re=-(w).c2.im; \
        (v).c2.im=(w).c2.re; \
        (v).c3.re=-(w).c3.im; \
        (v).c3.im=(w).c3.re;}

/* v=(-i)*w */
#define su3_vec_mul_mi(v,w){ \
        (v).c1.re=(w).c1.im; \
        (v).c1.im=-(w).c1.re; \
        (v).c2.re=(w).c2.im; \
        (v).c2.im=-(w).c2.re; \
        (v).c3.re=(w).c3.im; \
        (v).c3.im=-(w).c3.re;}

/* v=(-1)*w */
#define su3_vec_mul_m1(v,w){ \
        (v).c1.re=-(w).c1.re; \
        (v).c1.im=-(w).c1.im; \
        (v).c2.re=-(w).c2.re; \
        (v).c2.im=-(w).c2.im; \
        (v).c3.re=-(w).c3.re; \
        (v).c3.im=-(w).c3.im;}

/* r=s1+s2 */
#define su3vec_add(r,s1,s2){ \
        (r).c1.re=(s1).c1.re+(s2).c1.re; \
        (r).c1.im=(s1).c1.im+(s2).c1.im; \
        (r).c2.re=(s1).c2.re+(s2).c2.re; \
        (r).c2.im=(s1).c2.im+(s2).c2.im; \
        (r).c3.re=(s1).c3.re+(s2).c3.re; \
        (r).c3.im=(s1).c3.im+(s2).c3.im;}

/* r+=s */
#define su3vec_add_single(r,s){ \
        (r).c1.re+=(s).c1.re; \
        (r).c1.im+=(s).c1.im; \
        (r).c2.re+=(s).c2.re; \
        (r).c2.im+=(s).c2.im; \
        (r).c3.re+=(s).c3.re; \
        (r).c3.im+=(s).c3.im;}

/* r=s1-s2 */
#define su3vec_sub(r,s1,s2){ \
        (r).c1.re=(s1).c1.re-(s2).c1.re; \
        (r).c1.im=(s1).c1.im-(s2).c1.im; \
        (r).c2.re=(s1).c2.re-(s2).c2.re; \
        (r).c2.im=(s1).c2.im-(s2).c2.im; \
        (r).c3.re=(s1).c3.re-(s2).c3.re; \
        (r).c3.im=(s1).c3.im-(s2).c3.im;}

/* r-=s */
#define su3vec_sub_single(r,s){ \
        (r).c1.re-=(s).c1.re; \
        (r).c1.im-=(s).c1.im; \
        (r).c2.re-=(s).c2.re; \
        (r).c2.im-=(s).c2.im; \
        (r).c3.re-=(s).c3.re; \
        (r).c3.im-=(s).c3.im;}

/* r=c*s (c real) */
#define su3vec_real_mult(r,c,s){ \
        (r).c1.re=(c)*(s).c1.re; \
        (r).c1.im=(c)*(s).c1.im; \
        (r).c2.re=(c)*(s).c2.re; \
        (r).c2.im=(c)*(s).c2.im; \
        (r).c3.re=(c)*(s).c3.re; \
        (r).c3.im=(c)*(s).c3.im;}

/* r=u*s (u: SU(3) matrix; s: SU(3) vector) */
#define su3vec_su3_mult(r,u,s){ \
        (r).c1.re= (u).c11.re*(s).c1.re-(u).c11.im*(s).c1.im  \
                  +(u).c12.re*(s).c2.re-(u).c12.im*(s).c2.im  \
                  +(u).c13.re*(s).c3.re-(u).c13.im*(s).c3.im; \
        (r).c1.im= (u).c11.re*(s).c1.im+(u).c11.im*(s).c1.re  \
                  +(u).c12.re*(s).c2.im+(u).c12.im*(s).c2.re  \
                  +(u).c13.re*(s).c3.im+(u).c13.im*(s).c3.re; \
        (r).c2.re= (u).c21.re*(s).c1.re-(u).c21.im*(s).c1.im  \
                  +(u).c22.re*(s).c2.re-(u).c22.im*(s).c2.im  \
                  +(u).c23.re*(s).c3.re-(u).c23.im*(s).c3.im; \
        (r).c2.im= (u).c21.re*(s).c1.im+(u).c21.im*(s).c1.re  \
                  +(u).c22.re*(s).c2.im+(u).c22.im*(s).c2.re  \
                  +(u).c23.re*(s).c3.im+(u).c23.im*(s).c3.re; \
        (r).c3.re= (u).c31.re*(s).c1.re-(u).c31.im*(s).c1.im  \
                  +(u).c32.re*(s).c2.re-(u).c32.im*(s).c2.im  \
                  +(u).c33.re*(s).c3.re-(u).c33.im*(s).c3.im; \
        (r).c3.im= (u).c31.re*(s).c1.im+(u).c31.im*(s).c1.re  \
                  +(u).c32.re*(s).c2.im+(u).c32.im*(s).c2.re  \
                  +(u).c33.re*(s).c3.im+(u).c33.im*(s).c3.re;}

/* r=u^+*s (u: SU(3) matrix; s: SU(3) vector) */
#define su3vec_su3_dag_mult(r,u,s){ \
        (r).c1.re= (u).c11.re*(s).c1.re+(u).c11.im*(s).c1.im  \
                  +(u).c21.re*(s).c2.re+(u).c21.im*(s).c2.im  \
                  +(u).c31.re*(s).c3.re+(u).c31.im*(s).c3.im; \
        (r).c1.im= (u).c11.re*(s).c1.im-(u).c11.im*(s).c1.re  \
                  +(u).c21.re*(s).c2.im-(u).c21.im*(s).c2.re  \
                  +(u).c31.re*(s).c3.im-(u).c31.im*(s).c3.re; \
        (r).c2.re= (u).c12.re*(s).c1.re+(u).c12.im*(s).c1.im  \
                  +(u).c22.re*(s).c2.re+(u).c22.im*(s).c2.im  \
                  +(u).c32.re*(s).c3.re+(u).c32.im*(s).c3.im; \
        (r).c2.im= (u).c12.re*(s).c1.im-(u).c12.im*(s).c1.re  \
                  +(u).c22.re*(s).c2.im-(u).c22.im*(s).c2.re  \
                  +(u).c32.re*(s).c3.im-(u).c32.im*(s).c3.re; \
        (r).c3.re= (u).c13.re*(s).c1.re+(u).c13.im*(s).c1.im  \
                  +(u).c23.re*(s).c2.re+(u).c23.im*(s).c2.im  \
                  +(u).c33.re*(s).c3.re+(u).c33.im*(s).c3.im; \
        (r).c3.im= (u).c13.re*(s).c1.im-(u).c13.im*(s).c1.re  \
                  +(u).c23.re*(s).c2.im-(u).c23.im*(s).c2.re  \
                  +(u).c33.re*(s).c3.im-(u).c33.im*(s).c3.re;}



/* SU(3) fermion macros */

/* r=s1+s2 */
#define su3wferm_add(r,s1,s2){ \
        su3vec_add((r).d1,(s1).d1,(s2).d1); \
        su3vec_add((r).d2,(s1).d2,(s2).d2); \
        su3vec_add((r).d3,(s1).d3,(s2).d3); \
        su3vec_add((r).d4,(s1).d4,(s2).d4);}

/* r+=s */
#define su3wferm_add_single(r,s){ \
        su3vec_add_single((r).d1,(s).d1); \
        su3vec_add_single((r).d2,(s).d2); \
        su3vec_add_single((r).d3,(s).d3); \
        su3vec_add_single((r).d4,(s).d4);}

/* r=s1-s2 */
#define su3wferm_sub(r,s1,s2){ \
        su3vec_sub((r).d1,(s1).d1,(s2).d1); \
        su3vec_sub((r).d2,(s1).d2,(s2).d2); \
        su3vec_sub((r).d3,(s1).d3,(s2).d3); \
        su3vec_sub((r).d4,(s1).d4,(s2).d4);}

/* r-=s */
#define su3wferm_sub_single(r,s){ \
        su3vec_sub_single((r).d1,(s).d1); \
        su3vec_sub_single((r).d2,(s).d2); \
        su3vec_sub_single((r).d3,(s).d3); \
        su3vec_sub_single((r).d4,(s).d4);}

/* r=c*s (c real) */
#define su3wferm_real_mult(r,c,s){ \
        su3vec_real_mult((r).d1,c,(s).d1); \
        su3vec_real_mult((r).d2,c,(s).d2); \
        su3vec_real_mult((r).d3,c,(s).d3); \
        su3vec_real_mult((r).d4,c,(s).d4);}

/* v=gamma_0*w */
#define mul_su3wferm_g0(v,w){ \
        su3_vec_mul_m1((v).d1,(w).d3); \
        su3_vec_mul_m1((v).d2,(w).d4); \
        su3_vec_mul_m1((v).d3,(w).d1); \
        su3_vec_mul_m1((v).d4,(w).d2);}

/* v=gamma_1*w */
#define mul_su3wferm_g1(v,w){ \
        su3_vec_mul_mi((v).d1,(w).d4); \
        su3_vec_mul_mi((v).d2,(w).d3); \
        su3_vec_mul_i((v).d3,(w).d2); \
        su3_vec_mul_i((v).d4,(w).d1);}

/* v=gamma_2*w */
#define mul_su3wferm_g2(v,w){ \
        su3_vec_mul_m1((v).d1,(w).d4); \
        (v).d2=(w).d3; \
        (v).d3=(w).d2; \
        su3_vec_mul_m1((v).d4,(w).d1);}

/* v=gamma_3*w */
#define mul_su3wferm_g3(v,w){ \
        su3_vec_mul_mi((v).d1,(w).d3); \
        su3_vec_mul_i((v).d2,(w).d4); \
        su3_vec_mul_i((v).d3,(w).d1); \
        su3_vec_mul_mi((v).d4,(w).d2);}

/* v=-gamma_0*w */
#define mul_su3wferm_mg0(v,w){ \
        (v).d1=(w).d3; \
        (v).d2=(w).d4; \
        (v).d3=(w).d1; \
        (v).d4=(w).d2;}

/* v=-gamma_1*w */
#define mul_su3wferm_mg1(v,w){ \
        su3_vec_mul_i((v).d1,(w).d4); \
        su3_vec_mul_i((v).d2,(w).d3); \
        su3_vec_mul_mi((v).d3,(w).d2); \
        su3_vec_mul_mi((v).d4,(w).d1);}

/* v=-gamma_2*w */
#define mul_su3wferm_mg2(v,w){ \
        (v).d1=(w).d4; \
        su3_vec_mul_m1((v).d2,(w).d3); \
        su3_vec_mul_m1((v).d3,(w).d2); \
        (v).d4=(w).d1;}

/* v=-gamma_3*w */
#define mul_su3wferm_mg3(v,w){ \
        su3_vec_mul_i((v).d1,(w).d3); \
        su3_vec_mul_mi((v).d2,(w).d4); \
        su3_vec_mul_mi((v).d3,(w).d1); \
        su3_vec_mul_i((v).d4,(w).d2);}

/* v=gamma_5*w */
#define mul_su3wferm_g5(v,w){ \
        (v).d1=(w).d1; \
        (v).d2=(w).d2; \
        su3_vec_mul_m1((v).d3,(w).d3); \
        su3_vec_mul_m1((v).d4,(w).d4);}

/* v=u*w (u: SU(2) matrix; w: fermion) */
#define su3wferm_su3_mult(v,u,w){ \
        su3vec_su3_mult((v).d1,u,(w).d1); \
        su3vec_su3_mult((v).d2,u,(w).d2); \
        su3vec_su3_mult((v).d3,u,(w).d3); \
        su3vec_su3_mult((v).d4,u,(w).d4);}

/* v=u^+*w (u: SU(2) matrix; w: fermion) */
#define su3wferm_su3_dag_mult(v,u,w){ \
        su3vec_su3_dag_mult((v).d1,u,(w).d1); \
        su3vec_su3_dag_mult((v).d2,u,(w).d2); \
        su3vec_su3_dag_mult((v).d3,u,(w).d3); \
        su3vec_su3_dag_mult((v).d4,u,(w).d4);}

