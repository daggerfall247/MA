
/*******************************************************************************
*
* File misc.h
*
* Copyright (C) 2016 Bastian Brandt
* [Copyright (C) 2005 Martin Luescher -- source for some commands]
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Contains necessary definitions for the use of the routines from endian.c.
*
*******************************************************************************/

#define MISC_H

#include <limits.h>
#include <float.h>

#if ((DBL_MANT_DIG!=53)||(DBL_MIN_EXP!=-1021)||(DBL_MAX_EXP!=1024))
#error : Machine is not compliant with the IEEE-754 standard
#endif

#if (SHRT_MAX==0x7fffffff)
typedef short int stdint_t;
typedef unsigned short int stduint_t;
#elif (INT_MAX==0x7fffffff)
typedef int stdint_t;
typedef unsigned int stduint_t;
#elif (LONG_MAX==0x7fffffff)
typedef long int stdint_t;
typedef unsigned long int stduint_t;
#else
#error : There is no four-byte integer type on this machine 
#endif

#undef UNKNOWN_ENDIAN
#undef LITTLE_ENDIAN
#undef BIG_ENDIAN

#define UNKNOWN_ENDIAN 0
#define LITTLE_ENDIAN 1
#define BIG_ENDIAN 2

#ifndef ENDIAN_C
extern int endianness(void);
extern void bswap_int(int n,stdint_t *a);
extern void bswap_double(int n,double *a);
#endif
