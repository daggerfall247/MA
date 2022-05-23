
/*******************************************************************************
*
* File IO_utils.c
*
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes basic IO routines.
* 
* Externally accessile functions:
* 
* void logging(char *format,...)
*      Writes any type of information into the logfile.
*      Note: The argument consists of a format string and possible additional
*            entries.
*            -> usage like printf.
* 
* void error(int test,char *name,char *format,...)
*      Tests whether test=TRUE. If so, it aborts the program after printing
*      'name' and the message defines by the format string 'format' to the
*      log file.
*      Note: format and the following arguments work as for logging!
*            name should conventionally contain the name of the calling
*                 function and in brackets the filename where it is included.
*            -> Makes error detection easier!
* 
* void checkpoint(char *format,...)
*      In DEBUG-mode it writes a message to the logfile (like 'logging').
*      This is mostly used to indicate which function one is currently
*      entering.
*      -> Helps for bug detection!
*
*******************************************************************************/

#define IO_UTILS_C

#include<stdio.h>
#include<stdarg.h>
#include<stdlib.h>
#include"headers.h"
#include"modules.h"


void logging(char *format,...)
{
   FILE *flog=NULL;
   va_list args;

   flog=fopen(LOG_FILE,"ab");
   va_start(args,format);
   vfprintf(flog,format,args);
   va_end(args);
   fflush(flog);
   fclose(flog);
}


void error(int test,char *name,char *format,...)
{
   FILE *flog=NULL;
   va_list args;

   if(test!=0)
   {
      flog=fopen(LOG_FILE,"ab");
      if(flog==NULL)
         flog=fopen(LOG_FILE,"w");
      fprintf(flog,"\nError in %s:\n",name);
      va_start(args,format);
      vfprintf(flog,format,args);
      va_end(args);
      fprintf(flog,"\nProgram aborted!!!\n\n");
      fflush(flog);
      fclose(flog);
      exit(-1);
   }
}


void checkpoint(char *format,...)
{
#if(DEBUG==1)
   FILE *flog=NULL;
   va_list args;

   flog=fopen(LOG_FILE,"ab");
   if(flog==NULL)
      flog=fopen(LOG_FILE,"w");
   error(flog==NULL,0,"checkpoint [IO_utils.c]",
         "Cannot open logfile!");
   va_start(args,format);
   fprintf(flog,"\nEntering:\t");
   vfprintf(flog,format,args);
   va_end(args);
   fprintf(flog,"\n");
   fclose(flog);
#endif
}

