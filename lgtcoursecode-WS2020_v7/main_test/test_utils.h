
/*******************************************************************************
*
* File test.c
*
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2021 Alessandro Sciarra
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Utility functions for testing
*
*******************************************************************************/

#define MAIN_C

#include"ranlxd.h"
#include"modules.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

#define TEST_DEBUG 0
const double doubleRelativePrecision=1.e-12;
unsigned int testCounter=0;
enum fermionField {f_pointSource, f_cold, f_ascending};
enum gaugeField {g_cold, g_ascending};
enum resultType {sum, squareNorm};

typedef struct {
    enum fermionField fF;
    enum gaugeField gF;
    double mW;
    enum resultType rT;
} testParameters;


void setSpinorsToCold(sun_wferm *f)
{
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    for(;s<sf;s+=2)
    {
       *s=1.;
       *(s+1)=0.;
    }
}

void setSpinorsToAscending(sun_wferm *f)
{
    int ii;
    double *s,*sf;
    s=(double*)(f);
    sf=s+8*SUN*VOL;
    for(;s<sf;)
       for(ii=0;ii<8*SUN;ii++,s++)
          *s=(double)(ii+1);
}

void setGaugeFieldToCold(void)
{
   int n,dir;
   for(n=0;n<VOL;n++)
      for(dir=0;dir<DIM;dir++)
         sun_unit(*pu[n][dir]);
}

void setGaugeFieldToAscending(void)
{
   int n,dir,ii;
   double *du;
   for(n=0;n<VOL;n++)
      for(dir=0;dir<DIM;dir++)
      {
         du=(double*)(pu[n][dir]);
         for(ii=0;ii<SUNVOL;ii++,du++)
            *du=(double)(ii+1);
      }
}

bool isEqual(const double lhs, const double rhs, const double epsilon)
{
    if (lhs == 0.0 || rhs == 0.0)
        return fabs(lhs-rhs)<epsilon;
    else
        return fabs(lhs-rhs)/fabs(lhs)<=epsilon  &&  fabs(lhs-rhs)/fabs(rhs)<=epsilon;
}

bool makeTest_real(const double result, const double refValue, const char testType, const double precisionInPercent)
{
#if(TEST_DEBUG!=0)
    printf("   %25.15f == %25.15f   ", result, refValue);
#endif
    if (testType=='e')
    {
        double precision = precisionInPercent;
        if (precisionInPercent == 0.0)
            precision = doubleRelativePrecision;
        return isEqual(result, refValue, precision);
    }
    else
        return false;
}

bool makeTest_complex(const complex result, const complex refValue, const char testType, const double precisionInPercent)
{
    if (testType=='e')
    {
        bool real = makeTest_real(result.re, refValue.re, 'e', precisionInPercent);
        bool imag = makeTest_real(result.im, refValue.im, 'e', precisionInPercent);
        return  real && imag;
    }
    else
        return false;
}

void printHeadlines(void)
{
    printf("\n ---------------------------------------------------------------------------------\n");
    printf("   Running tests...\n");
}

void printReportLine(const bool passed, const char* testName, const unsigned int testNumber)
{
    printf("     %3d  %-60s   %s\n", testNumber, testName, (passed) ? "\033[92mOK\033[0m" : "\033[91mFAILED\033[0m");
}

void printFootlines(const unsigned int passed, const unsigned int failed)
{
    printf("   ...done!");
    printf("\n ---------------------------------------------------------------------------------\n");
    printf("        \033[92m%d test(s) passed\033[0m, \033[91m%d test(s) failed\033[0m!", passed, failed);
    printf("\n ---------------------------------------------------------------------------------\n\n");
}

double getReferenceValue_applyWilsonDiracOperator(testParameters tP)
{
    if (tP.mW == 1.0)
    {
        if (tP.fF == f_cold && tP.gF == g_cold)
        {
            if (tP.rT == squareNorm)
                return 3.0;
            else if (tP.rT == sum)
                return 1.5;
        }
        else if (tP.fF == f_ascending && tP.gF == g_cold)
        {
            if (tP.rT == squareNorm)
                return 1225.0;
            else if (tP.rT == sum)
                return 37.5;
        }
        else if (tP.fF == f_cold && tP.gF == g_ascending)
        {
            if (tP.rT == squareNorm)
                return 14677.0;
            else if (tP.rT == sum)
                return -104.5;
        }
        else if (tP.fF == f_ascending && tP.gF == g_ascending)
        {
            if (tP.rT == squareNorm)
                return 6098236.33333333;
            else if (tP.rT == sum)
                return -1966.5;
        }
    }
#if(TEST_DEBUG!=0)
    else
        fprintf(stderr,"ERROR: No reference value known for given test!\n");
#endif

    return nan("");
}

double calculateTestResult(sun_wferm *f, testParameters tP)
{
    if (tP.rT == sum)
        return globalSum(f);
    else if (tP.rT == squareNorm)
        return globalSquareNorm(f);
#if(TEST_DEBUG!=0)
    else
        fprintf(stderr,"ERROR: No procedure known to calculate result of given test!\n");
#endif

    return nan("");
}

void initializeFermionFieldBasedOnTestParameters(testParameters tP, sun_wferm *f1)
{
    if (tP.fF == f_pointSource)
        pointSource(f1, 0, 0);
    else if (tP.fF == f_cold)
        setSpinorsToCold(f1);
    else if (tP.fF == f_ascending)
        setSpinorsToAscending(f1);
#if(TEST_DEBUG!=0)
    else
    {
        fprintf(stderr,"ERROR: Invalid parameters for to testApplyWilsonDiracOperator!");
        f1 = NULL;
    }
#endif
}

void initializeGaugeFieldBasedOnTestParameters(testParameters tP)
{
    if (tP.gF == g_cold)
        setGaugeFieldToCold();
    else if (tP.gF == g_ascending)
        setGaugeFieldToAscending();
#if(TEST_DEBUG!=0)
    else
    {
        fprintf(stderr,"ERROR: Invalid parameters for to testApplyWilsonDiracOperator!");
        pu[0][0] = NULL;
    }
#endif
}

bool testApplyWilsonDiracOperator(testParameters tP)
{
    runParams.mass=tP.mW;

    sun_wferm *f1,*f2;
    allocateFermionField(&f1);
    allocateFermionField(&f2);

    initializeFermionFieldBasedOnTestParameters(tP, f1);
    initializeGaugeFieldBasedOnTestParameters(tP);
#if(TEST_DEBUG!=0)
    if (pu[0][0]==NULL || f1==NULL)
    {
        deallocateFermionField(&f1);
        deallocateFermionField(&f2);
        return false;
    }
#endif

    applyWilsonDiracOperator(f2,f1);
    bool toBeReturned = makeTest_real(getReferenceValue_applyWilsonDiracOperator(tP), calculateTestResult(f2, tP), 'e', 0);
    deallocateFermionField(&f1);
    deallocateFermionField(&f2);
    return toBeReturned;
}

bool testApplyWilsonDiracOperatorDagger(testParameters tP)
{
    runParams.mass=tP.mW;

    complex DpsidagDpsi, psidagDdagDpsi;
    sun_wferm *s, *f1,*f2;
    allocateFermionField(&s);
    allocateFermionField(&f1);
    allocateFermionField(&f2);

    initializeFermionFieldBasedOnTestParameters(tP, f1);
    initializeGaugeFieldBasedOnTestParameters(tP);
#if(TEST_DEBUG!=0)
    if (pu[0][0]==NULL || f1==NULL)
    {
        deallocateFermionField(&s);
        deallocateFermionField(&f1);
        deallocateFermionField(&f2);
        return false;
    }
#endif

    applyWilsonDiracOperator(f2,f1);
    DpsidagDpsi = scalarProd(f2,f2);
    applyComplexConjDiracOperator(s,f2);
    psidagDdagDpsi = scalarProd(f1,s);

    bool toBeReturned = makeTest_complex(DpsidagDpsi, psidagDdagDpsi, 'e', 0);
    deallocateFermionField(&s);
    deallocateFermionField(&f1);
    deallocateFermionField(&f2);
    return toBeReturned;
}

bool testCGSolver(testParameters tP)
{
    runParams.mass=tP.mW;

    sun_wferm *s, *f1,*f2;
    allocateFermionField(&s);
    allocateFermionField(&f1);
    allocateFermionField(&f2);

    initializeFermionFieldBasedOnTestParameters(tP, s);
    initializeGaugeFieldBasedOnTestParameters(tP);
#if(TEST_DEBUG!=0)
    if (pu==NULL || s==NULL)
        return false;
#endif

    applyComplexConjDiracOperator(f2,s);
    allocateFermionFieldsForCG();
    setAllSpinorsToZero(f1);
    int cgSteps=cg(f1,applyWilsonDiracOperator,applyComplexConjDiracOperator,f2,1.e-11,1000);
    if (cgSteps<0)
    {
#if(TEST_DEBUG!=0)
        fprintf(stderr,"ERROR: CG did not converge in testCGSolver!");
#endif
        return false;
    }

    bool chi_test  = makeTest_real(4*SUN*VOL*realPartOfScalarProd(s,f1), 1.659500e-01, 'e', 1.e-6);
    applyWilsonDiracOperator(f2,f1);
    multiplyByRealAndSum(f1,f2,-1.,s);
    bool norm_test = makeTest_real(sqrt(4*SUN*VOL*globalSquareNorm(f1)), 0.0, 'e', 0);
    deallocateFermionFieldsForCG();
    deallocateFermionField(&s);
    deallocateFermionField(&f1);
    deallocateFermionField(&f2);
    return chi_test && norm_test;
}

void getReferenceValue_test2ptFunction(testParameters tP, double *refValue)
{
    if (tP.mW == 2.0 && tP.fF == f_pointSource && tP.gF == g_cold)
    {
        refValue[0] = 3.64692520e-01;
        refValue[1] = 7.05856425e-03;
        refValue[2] = 6.93918804e-04;
        refValue[3] = 7.05856425e-03;
        return;
    }
#if(TEST_DEBUG!=0)
    else
        fprintf(stderr,"ERROR: No reference value known for given test!\n");
#endif

    for(int i=0; i<LENGT; i++)
        refValue[i]=nan("");
}

void parseLogFileToGet2ptFunctionValues(double *result)
{
    FILE *flog=fopen(LOG_FILE,"r");
    if(flog==NULL)
        error(1,"test_2pt [test_utils.c]","Unable to open logfile!");
    char * line = NULL;
    size_t bufsize = 32;
    line = (char *)malloc(bufsize * sizeof(char));
    if( line == NULL)
    {
        perror("Unable to allocate buffer");
        exit(1);
    }
    size_t len = 0;
    ssize_t read;
    int dataRead = 0, timeslice;
    double value;
    while ((read = getline(&line, &len, flog)) != -1 && dataRead<4) {
        if(strncmp(line, "* 2pt :", 7) == 0){
            sscanf (line,"* 2pt : %*d %*d %*d %d %lf %*f", &timeslice, &value);
            result[timeslice]=value;
            dataRead++;
        }
    }
    remove(LOG_FILE);

}

bool test2ptFunction(testParameters tP)
{
    int gammaToBeUsed=5;
    char out_dir[] = ".";
    runParams.mass=tP.mW;
    measParams.source[0]=0;
    measParams.source[1]=0;
    measParams.source[2]=0;
    measParams.source[3]=0;
    measParams.stype=1;
    measParams.nmax=5000;
    measParams.eps=1.e-8;
    measParams.n2pt=1;
    measParams.g1=&gammaToBeUsed;
    measParams.g2=&gammaToBeUsed;
    setupOutputFilesForMeas(runParams.idForOutputFilesName,out_dir);

    allocateFermionFieldsFor2ptFunctions(measParams.stype);
    initializeGaugeFieldBasedOnTestParameters(tP);
#if(TEST_DEBUG!=0)
    if (pu==NULL)
        return false;
#endif

    measure2ptFunctions(measParams.source,measParams.stype,measParams.n2pt,measParams.g1,measParams.g2);

    double correlatorValues[4];
    parseLogFileToGet2ptFunctionValues(correlatorValues);

    double reference[4];
    getReferenceValue_test2ptFunction(tP, reference);

    bool testResult = true;
    for(int i=0; i<4; i++)
        testResult = testResult && makeTest_real(correlatorValues[i], reference[i], 'e', 0);

    deallocateFermionFieldsFor2ptFunctions();
    remove(LOG_FILE);
    return testResult;
}

bool runSingleTest(bool (*testFunction)(testParameters), testParameters tP, const char* testName){
    bool outcome=testFunction(tP);
    testCounter++;
    printReportLine(outcome, testName, testCounter);
    return outcome;
}

void runAllTests(void){
    unsigned int passed=0;
    printHeadlines();

    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_cold,      g_cold,      1.0, squareNorm},    "D.psi [f_cold,      g_cold,     1.0, squareNorm]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_cold,      g_cold,      1.0, sum       },    "D.psi [f_cold,      g_cold,     1.0, sum       ]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_ascending, g_cold,      1.0, squareNorm},    "D.psi [f_ascending, g_cold,     1.0, squareNorm]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_ascending, g_cold,      1.0, sum       },    "D.psi [f_ascending, g_cold,     1.0, sum       ]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_cold,      g_ascending, 1.0, squareNorm},    "D.psi [f_cold,      g_ascending 1.0, squareNorm]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_cold,      g_ascending, 1.0, sum       },    "D.psi [f_cold,      g_ascending 1.0, sum       ]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_ascending, g_ascending, 1.0, squareNorm},    "D.psi [f_ascending, g_ascending 1.0, squareNorm]");
    passed+=runSingleTest(testApplyWilsonDiracOperator, (testParameters){f_ascending, g_ascending, 1.0, sum       },    "D.psi [f_ascending, g_ascending 1.0, sum       ]");

    passed+=runSingleTest(testApplyWilsonDiracOperatorDagger, (testParameters){f_cold,      g_cold,      1.0, squareNorm},    "Ddagger_test [f_cold,      g_cold,      1.0]");
    passed+=runSingleTest(testApplyWilsonDiracOperatorDagger, (testParameters){f_ascending, g_cold,      1.0, squareNorm},    "Ddagger_test [f_ascending, g_cold,      1.0]");
    passed+=runSingleTest(testApplyWilsonDiracOperatorDagger, (testParameters){f_cold,      g_ascending, 1.0, squareNorm},    "Ddagger_test [f_cold,      g_ascending, 1.0]");
    passed+=runSingleTest(testApplyWilsonDiracOperatorDagger, (testParameters){f_ascending, g_ascending, 1.0, squareNorm},    "Ddagger_test [f_ascending, g_ascending, 1.0]");

    passed+=runSingleTest(testCGSolver, (testParameters){f_pointSource, g_cold, 2.0, squareNorm},    "chiral condensate [f_pointSource, g_cold, 2.0]");

    passed+=runSingleTest(test2ptFunction, (testParameters){f_pointSource, g_cold, 2.0, squareNorm},    "2-point function [f_pointSource, g_cold, 2.0]");

    printFootlines(passed, testCounter-passed);
}


