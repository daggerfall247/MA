
/*******************************************************************************
*
* File modules.h
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Provides the prototypes for the global functions from all files.
*
*******************************************************************************/


#ifndef GAUGE_H
 #include"gauge.h"
#endif

#ifndef HEADERS_H
 #include"headers.h"
#endif

/* Observables */

#ifndef PLAQUETTE_C
extern double plaquette(void);
extern double gaugeAction(void);
#endif

#ifndef SMEARING_C
extern void smearing_APE_spatial(int,double,sun_mat *u[VOL][DIM]);
extern void smearing_APE_temporal(int,double,sun_mat *u[VOL][DIM]);
#endif

#ifndef WILSON_C
extern void measureWilsonLoop(double*);
#endif

#ifndef TPT_C
extern void allocateFermionFieldsFor2ptFunctions(int);
extern void deallocateFermionFieldsFor2ptFunctions(void);
extern void measure2ptFunctions(int*,int,int,int*,int*);
#endif

/* Dirac operators */
#ifndef DIRAC_WIL_C
extern void applyWilsonDiracOperator(sun_wferm*,sun_wferm*);
extern void applyComplexConjDiracOperator(sun_wferm*,sun_wferm*);
#endif

#ifndef SPIN_ALG_C
extern void setAllSpinorsToZero(sun_wferm*);
extern void copySpinors(sun_wferm*,sun_wferm*);
extern void multiplyByRealAndSum(sun_wferm*,sun_wferm*,double,sun_wferm*);
extern double globalSquareNorm(sun_wferm*);
extern double globalSum(sun_wferm*);
extern double realPartOfScalarProd(sun_wferm*,sun_wferm*);
extern complex scalarProd(sun_wferm *f2,sun_wferm *f1);
#endif

#ifndef SOURCES_C
extern void pointSource(sun_wferm*,int,int);
#endif

/* Solver */

#ifndef SOLV_CG_C
extern void allocateFermionFieldsForCG(void);
extern void deallocateFermionFieldsForCG(void);
extern int cg(sun_wferm*,void (*A)(sun_wferm *r,sun_wferm *s),
              void (*Ad)(sun_wferm *r,sun_wferm *s),sun_wferm*,
              double,int);
#endif

/* Updates */

#ifndef EXP_FCT_C
extern void expx(double,sun_alg*,sun_mat*);
#endif

#ifndef METRO_C
extern void staples(int,int,sun_mat*);
extern double localMetropolisUpdate(int,int,int);
#endif

#ifndef HEATBATH_C
extern void staples(int,int,sun_mat*);
extern double localHeatbathUpdate(int,int,int);
#endif

#ifndef UPDATE_C
extern void gaugefieldUpdate(int,int,int,int);
#endif

/* IO-archive functions */

#ifndef INP_IO_C
extern void readInputFile(int*,char*,int*,char*,int,char *argv[]);
extern void readInputFileForMeas(int *,int *,int *,int *,char *,int,char *argv[]);
extern void setupOutputFiles(int,char*);
extern void setupOutputFilesForMeas(int,char*);
extern void printStartupInfo(int,int,char*);
#endif

#ifndef IO_UTILS_C
extern void logging(char*,...);
extern void error(int,char*,char*,...);
extern void checkpoint(char *format,...);
#endif

#ifndef CONFIG_IO_C
extern void writeConfig(char*);
extern void readConfig(char*);
#endif

/* Initialisation */

#ifndef INIT_C
extern void initProgram(int);
extern void initArrayOfNeighbours(void);
extern void initGaugeField(int);
extern void releaseGaugeField(void);
extern void allocateGaugeField(sun_mat *u[VOL][DIM]);
extern void deallocateGaugeField(sun_mat *u[VOL][DIM]);
extern void copyGaugeField(sun_mat *u1[VOL][DIM],sun_mat *u2[VOL][DIM]);
extern void allocateFermionField(sun_wferm**);
extern void deallocateFermionField(sun_wferm**);
#endif

#ifndef RANDOM_SU3_C
extern void initTwoPi(void);
extern int checkTwoPi(void);
extern void gauss(double*,int);
extern void su3RandomVector(su3vec*);
extern void su3RandomMatrix(su3mat*);
extern void su2RandomMatrix(su2mat*);
#endif

/* utils functions */

#ifndef UTILS_C
extern double getTime(void);
extern int custom_isnan(double);
extern int custom_isinf(double);
#endif

#ifndef ERROR_CHECKS_C
extern void setWarning(void);
extern void checkForErrors(int,int);
#endif

/* maths */

#ifndef SUN_VFUNC_C
extern void project_to_su3(su3mat *u);
extern void project_to_su2(su2mat *u);
extern void project_gfield_to_sun(sun_mat *u[VOL][DIM]);
#endif
