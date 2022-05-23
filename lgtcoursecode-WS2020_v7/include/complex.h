
/*******************************************************************************
*
* File complex.h
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Includes definitions and macros for complex numbers
*
*******************************************************************************/

#define COMPLEX_H

#ifndef MATH_H
 #include<math.h>
#endif

typedef struct
{
   double re,im;
} complex;

/* Complex macros */

/* u=v*w */
#define compl_mult(u,v,w){ \
        (u).re=(v).re*(w).re-(v).im*(w).im; \
        (u).im=(v).re*(w).im+(v).im*(w).re;}

/* Re(u=v*w) */
#define compl_mult_re(u,v,w){ \
        (u)=(v).re*(w).re-(v).im*(w).im;}

/* n=|u| */
#define compl_norm(n,u){ \
        (n)=sqrt((u).re*(u).re+(u).im*(u).im);}

/* p=phase(u) */
#define get_phase(p,u){ \
        (p)=atan2((u).im,(u).re);}

/* u=star(v)*w */
#define compl_mult_sn(u,v,w){ \
        (u).re=(v).re*(w).re+(v).im*(w).im; \
        (u).im=(v).re*(w).im-(v).im*(w).re;}

/* Re(u=star(v)*w) */
#define compl_mult_sn_re(u,v,w){ \
        (u)=(v).re*(w).re+(v).im*(w).im;}

/* u=v+w */
#define compl_add(u,v,w){ \
        (u).re=(v).re+(w).re; \
        (u).im=(v).im+(w).im;}

/* u=v-w */
#define compl_sub(u,v,w){ \
        (u).re=(v).re-(w).re; \
        (u).im=(v).im-(w).im;}

/* u=u+v */
#define compl_selfadd(u,v){ \
        (u).re+=(v).re; \
        (u).im+=(v).im;}

/* u=u-v */
#define compl_selfsub(u,v){ \
        (u).re-=(v).re; \
        (u).im-=(v).im;}

/* u=u/v : v real */
#define compl_realdiv_single(u,v){ \
        (u).re/=(v); \
        (u).im/=(v);}

/* u=u*v : v real */
#define compl_realmul_single(u,v){ \
        (u).re*=(v); \
        (u).im*=(v);}
        
/* u=v/w */
#define compl_div(u,v,w){ \
        double zwii; \
        compl_mult_sn(u,w,v); \
        zwii=(w).re*(w).re+(w).im*(w).im; \
        (u).re/=zwii; \
        (u).im/=zwii;}

/* following as above with w=i*w */

/* u=star(v)*i*w */
#define compl_im_x_mult_sn(u,v,w){ \
        (u).re=-(v).re*(w).im+(v).im*(w).re; \
        (u).im=(v).re*(w).re+(v).im*(w).im;}

/* Re(u=star(v)*i*w) */
#define compl_im_x_mult_sn_re(u,v,w){ \
        (u)=-(v).re*(w).im+(v).im*(w).re;}

/* u=v+i*w */
#define compl_im_x_add(u,v,w){ \
        (u).re=(v).re-(w).im; \
        (u).im=(v).im+(w).re;}

/* u=v-i*w */
#define compl_im_x_sub(u,v,w){ \
        (u).re=(v).re+(w).im; \
        (u).im=(v).im-(w).re;}

/* u=u+i*w */
#define compl_im_x_selfadd(u,v){ \
        (u).re-=(v).im; \
        (u).im+=(v).re;}

/* u=u-i*w */
#define compl_im_x_selfsub(u,v){ \
        (u).re+=(v).im; \
        (u).im-=(v).re;}

/* u=star(v)+w */
#define compl_add_sn(u,v,w){ \
        (u).re=(v).re+(w).re; \
        (u).im=-(v).im+(w).im;}

/* u=star(v)-w */
#define compl_sub_sn(u,v,w){ \
        (u).re=(v).re-(w).re; \
        (u).im=-(v).im-(w).im;}
