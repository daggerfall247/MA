***************************************************************************
* 
* File README.txt
* 
* Copyright (C) 2020 Alessandro Sciarra
* Copyright (C) 2019 Francesca Cuteri
* Copyright (C) 2016 Bastian Brandt
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
* 
* Contains information about the program structure.
* 
* For more questions about the program please write to:
* -> sciarra@itp.uni-frankfurt.de
* ->  cuteri@itp.uni-frankfurt.de
* 
***************************************************************************

* Compiling and running the program *
**************************************

The program has one directory (main), which, upon compilation via make,
contains the executable program (qcd)

The main directory contains a Makefile to be used to generate the
executable (qcd). The Makefile can be extended to include additional
files and folders.

To execute the program you can go to the main directory and type:
./qcd -i test.in

The file following '-i' is the input file which is needed to hand over
the necessary parameters to the program. A sample input file 'test.in'
can be found in the main directory.

* Code structure: *
*******************

The code is organised in different directories. Those are:

Basic directories:
include   : Includes all .h files.
modules   : Includes the .c files for the functions used in the program.
main      : Includes the main program 'qcd.c' and the Makefile

Subdirectories in modules:
modules/admin  : Contains files which include adminitrative programs,
                 such as programs to initialise fields and other arrays.
                 It also contains files with general tools/functions
                 that do not fit in any of the other categories.
modules/io     : Contains the input and output functions of the program.
modules/random : Contains the functions to construct random objects,
                 such as random SU(3) matrices etc.
                 The random number generator 'ranlux' is also included.

-> NOTE: While, some of the above-mentioned functionalities are already
   there, some others are to be developed, in the course of the tutorial,
   within the corresponding subdirectories.

* Files included in the modules directories *
*********************************************

modules/admin/endian.c       : Routines for the determination of the
                               endiannes of the computer and byte-swaps.
modules/admin/error_checks.c : Routines to check for bugs and missing
                               initialisations of parts of the program.
modules/admin/init.c         : Routines for the startup of the program
                               and the initialisation of fields.
modules/admin/sun_vfunc.c    : Some utils for SU(2,3) matrices.
modules/admin/utils.c        : Some general utils.

modules/io/inp_IO.c   : Routines to read the input files and the basic
                        setup for output files.
modules/io/IO_utils.c : Some utils for messages and error handling.

modules/random/random_su3.c : Routines to generate random SU(3) matrices
                              and similar objects.
modules/random/ranlxd.c     : 'ranlux' random number generator.


For an overview over the available functions see the file 'modules.h'.
-> NOTE: For every new functions that has to be externally accessible,
   the corresponding prototype has to be added in 'modules.h'.
A description of the included externally accessible functions is also
given in the file headers of the individual files.


* Structure of .h files and global variables *
**********************************************

include/complex.h    : Contains the definition of complex numbers used
                       in the program. It also contains macros for
                       complex numbers.
include/gauge.h      : Contains the definition of structures related to
                       SU(2) and SU(3) matrices.
                       It also contains macros to work with the structs.
    -> NOTE: You should have a good look at the included macros!
include/headers.h    : Main .h file. Includes the global variables and
                       global definitions of macros. It also contains
                       sanity checks for the parameters that have to be
                       specified at compile time. It also contains the
                       redefinitions of macro names for the switch
                       between SU(2) and SU(3) gauge theory.
    -> NOTE: The way this is set up here is a particular choice how a
       program can handle simulations with a general gauge group SU(N).
    -> NOTE: Only N=2,3 is possible, but the structure can be extended.
include/misc.h       : Contains the necessary macros and prototypes for
                       the endiannes programs.
include/modules.h    : Contains the prototypes of all the globally
                       available functions.


* Lattice setup *
*****************

Fields will be represented by 1d arrays of the associated type of
variables/structs. A particular lattice point is then given by
a 1d index n. Given the extents LENGT (LT), LENGS1 (L1), LENGS2 (L2),
LENGS3 (L3) the index for the point (x0,x1,x2,x3) is given by:
 n = L1*L2*L3*x0 + L2*L3*x1 + L3*x2 + x3

To move on the lattice it is useful to define so-called 'neighbour'
arrays:
int neib[VOL][2*DIM]
The entry neib[n][i] should contain the index (n) of the neighbouring
point in +i direction while the entry neib[n][i+DIM] contains the
index of the neighbouring point in -i direction. These arrays need
to be initialised once in the beginning of the program and should be
globally accessible.


* General remarks *
*******************

All files must begin with a header stating the filename and its purpose.
It is also useful to explain in short the globally accessible functions.
As a guideline how to set this up you can use the headers of the included
files.

It is useful to define macros in the form of the filename at the
beginning of each file. This serves to manage the definition of global
prototypes and variables in different files (cf. include guards).

It is useful to extend the present structure when the program keeps
growing. Hints on how to do this will be given in the exercises.
