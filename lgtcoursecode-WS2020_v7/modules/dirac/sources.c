
/*******************************************************************************
*
* File sources.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routine needed to generate sources for inversion.
* 
* Externally accessible functions:
* 
* 
* 
*******************************************************************************/

#define SOURCES_C

#include"headers.h"
#include"modules.h"



void pointSource(sun_wferm *f,int n,int d)
{
    double *s,*sf;
    error(d>=4*SUN,"point_source [sources.c]",
          "Dirac index too large!");
    s =(double*)f;
    sf=s+8*SUN*VOL;
    for(;s<sf;s+=1)
    {
       *s=0.;
    }
    s=(double*)(f)+8*SUN*n;
    s+=2*d;
    *s=1.;
}
