
/*******************************************************************************
*
* File test.c

* Copyright (C) 2021 Alessandro Sciarra
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Program under construction for testing lattice gauge theory code
*
* Syntax: test -i <input-filename>
*
* For usage instructions see the file README.main
*
*******************************************************************************/

#define MAIN_C

#include"ranlxd.h"
#include"modules.h"
#include"test_utils.h"
#include <assert.h>

int main(int argc,char *argv[])
{
   static_assert (LENGT==4 || LENGS1==4 || LENGS2==4 || LENGS3==4,
           "To compile and run tests, please set lattice to 4^4.");
   initArrayOfNeighbours();
   initGaugeField(1);
   runAllTests();

   return 0;
}
