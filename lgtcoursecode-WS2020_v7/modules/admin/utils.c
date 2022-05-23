
/*******************************************************************************
*
* File utils.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes some additional routines.
* 
* Externally accessible functions:
* 
* double getTime(void)
*        Returns the current execution time after the first call of getTime.
*        I.e. the first call return 0!
* 
* int custom_isnan(double var)
*     Returns true when var is not-a-number.
* 
* int custom_isinf(double x)
*     Returns true if x is infinite.
* 
*******************************************************************************/

#define UTILS_C

#include <sys/time.h>
#include <stddef.h>

static int time_flag=0;
/*static time_t ref_time;*/
struct timeval t1;


//Unfortunately, the precision of this routine is just seconds. For a Nanosecond precision, see the commented parts. However, they might not work on all OS!
/*double getTime(void)*/
/*{*/
/*   if(time_flag==0)*/
/*   {*/
/*      time(&ref_time);*/
/*      time_flag=1;*/
/*      return 0.;*/
/*   }*/
/*   else*/
/*   {*/
/*      time_t wt;*/
/*      time(&wt);*/
/*      return difftime(wt,ref_time);*/
/*   }*/
/*}*/


double getTime(void)
{
   if(time_flag==0)
   {
      /*time(&ref_time);*/
      gettimeofday(&t1, 0);
      time_flag=1;
      return 0.;
   }
   else
   {
      /*time_t wt;
      time(&wt);
      return difftime(wt,ref_time);*/
      struct timeval t2;
      double dt;
      gettimeofday(&t2, 0);
      dt = t2.tv_sec - t1.tv_sec;
      dt += (t2.tv_usec - t1.tv_usec)/1.e6;
      return dt;
   }
}


int custom_isnan(double var)
{
    volatile double d = var;
    return d != d;
}


int custom_isinf(double x)
{
   volatile double temp = x;
   if ((temp == x) && ((temp - x) != 0.0))
      return 1;
   else return 0;
}
