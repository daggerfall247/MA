
/*******************************************************************************
*
* File 2pt.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Includes the routines to compute 2pt functions.
* 
*******************************************************************************/

#define TPT_C

#include<stdio.h>
#include<stdlib.h>
#include"modules.h"

static int init=0;
static int nsource;
static sun_wferm *s,*f1,*f2;
static sun_wferm **d;



static int get_index(int *x)
{
   int ii,n;
   for(ii=0,n=0;ii<DIM;ii++)
      n+=latParams.volumeFactor[ii]*x[ii];
   return n;
}


static void mk_chi_g5_point(int n,sun_wferm **d,complex **chi)
{
   int ia,ib;
   complex z,prod;
   complex *s1,*s2,*sf;

   for(ia=0;ia<nsource;ia++)
   {
      for(ib=0;ib<nsource;ib++)
      {
         s1=(complex*)&(d[ia][n]);
         s2=(complex*)&(d[ib][n]);
         sf=s1+4*SUN;
         prod.re=0.;
         prod.im=0.;
         for(;s1<sf;s1+=1,s2+=1)
         {
             compl_mult_sn(z,*s1,*s2);
             compl_selfadd(prod,z);
         }
         chi[ia][ib]=prod;
      }
   }
}


static complex contr_chi_g5(int nsource,complex **chi)
{
   int ia;
   complex sum;

   sum.re=0.;
   sum.im=0.;
   for(ia=0;ia<nsource;ia++)
      compl_selfadd(sum, chi[ia][ia]);

   return sum;
}


static void write_2pt(int n2pt,complex **c,int *g1,int *g2)
{
   FILE *flog=NULL;
   int ic,it;

   flog=fopen(LOG_FILE,"ab");
   if(flog==NULL)
      error(1,"write_2pt [2pt.c]","Unable to open logfile!");

   fprintf(flog,"\nMeasured 2pt functions:\n\n");
   for(ic=0;ic<n2pt;ic++)
   {
      for(it=0;it<latParams.linearExtent[0];it++)
      {
         fprintf(flog,"* 2pt : %d %d %d %d %.8e %.8e\n",
            ic,g1[ic],g2[ic],it,c[ic][it].re,c[ic][it].im);
      }
   }

   fflush(flog);
   fclose(flog);
}


void allocateFermionFieldsFor2ptFunctions(int stype)
{
   int ii;

   error(init!=0,"init_2pt [2pt.c]","Already initialised!");

   allocateFermionField(&s);
   allocateFermionField(&f1);
   allocateFermionField(&f2);

   if(stype==1)
      nsource=4*SUN;
   else
      error(1,"init_2pt [2pt.c]","Unknown source type!");
   d=malloc(nsource*sizeof(sun_wferm*));
   for(ii=0;ii<nsource;ii++)
      allocateFermionField(d+ii);

   allocateFermionFieldsForCG();
   init=1;
}


void deallocateFermionFieldsFor2ptFunctions(void)
{
   int ii;

   error(init==0,"finish_2pt [2pt.c]","Has not been initialised!");

   for(ii=0;ii<nsource;ii++)
      deallocateFermionField(d+ii);

   deallocateFermionField(&s);
   deallocateFermionField(&f1);
   deallocateFermionField(&f2);
}


void measure2ptFunctions(int *source,int stype,int n2pt,int *g1,int *g2)
{
   int ii,ic,n,test;
   int it,itt,nts,in,nn;
   double t1,t2,tt1,tt2;
   complex **ch,**c,z;
   n=get_index(source);
   error(init==0,"meas_2pt [2pt.c]","Has not been initialised!");

   t1=getTime();
   for(ii=0;ii<nsource;ii++)
   {
      if(stype==1)
         pointSource(s,n,ii);
      setAllSpinorsToZero(f1);
      applyComplexConjDiracOperator(f2,s);

      tt1=getTime();
      test=cg(f1,applyWilsonDiracOperator,applyComplexConjDiracOperator,f2,measParams.eps,measParams.nmax);
      tt2=getTime();
      error(test<0,"test","Inversion %d failed!",ii);
      logging("Inversion %d: %d iterations; took %.3f sec\n",ii,test,tt2-tt1);
      copySpinors(d[ii],f1);
   }

   ch=malloc(nsource*sizeof(complex*));
   for(ii=0;ii<nsource;ii++)
      ch[ii]=malloc(nsource*sizeof(complex));
   c=malloc(n2pt*sizeof(complex*));

   for(ic=0;ic<n2pt;ic++)
   {
      c[ic]=malloc(latParams.linearExtent[0]*sizeof(complex));
      for(it=0;it<latParams.linearExtent[0];it++)
      {
         itt=(it+source[0])%latParams.linearExtent[0];
         nts=itt*latParams.volumeOtherDirs[0];
         c[ic][it].re=0.;
         c[ic][it].im=0.;
         for(in=0;in<latParams.volumeOtherDirs[0];in++)
         {
            nn=nts+in;
            if(stype==1)
            {
               if(g1[ic]==5)
                  mk_chi_g5_point(nn,d,ch);
               else
                  error(1,"meas_2pt [2pt.c]","Unknown Gamma1!");
               if(g2[ic]==5)
                  z=contr_chi_g5(nsource,ch);
               else
                  error(1,"meas_2pt [2pt.c]","Unknown Gamma2!");
               compl_selfadd(c[ic][it],z);
            }
         }
      }
   }
   t2=getTime();

   logging("Computation of 2pt fct. took %.3f sec!\n",t2-t1);
   write_2pt(n2pt,c,g1,g2);

   for(ic=0;ic<n2pt;ic++)
      free(c[ic]);
   free(c);
   for(ii=0;ii<nsource;ii++)
      free(ch[ii]);
   free(ch);
}
