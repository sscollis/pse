# PSE:  Parabolized Stability Equations Solver

This solves the incompressible Parabolized Stability Equations (PSE).

While the code support both linear and nonlinear analysis, the nonlinear
capability is currently deactivated in the main `src` directory since it uses
either CRAY or SGI FFT routines that must now be replaced (ideally with FFTW).

There is also usage of deprecated LINPACK routines that should be replaced by
LAPACK.

Finally, the build uses a couple of commercially licensed routines from 
Numerical-Recipes in FORTRAN that are turned on using

## How to Build

    ln -s gcc.mak Makefile
    make USE_NR=1 

Note that the Numerical-Recipies code is not distributed with PSE and can only
be used by use under the terms of the commercial license.  It would be great
to strip out this NR code and it would be easy to do so...  This is left as
as project for a future developer.

## Directory structure

Directory    |  Description
-------------|------------------------------------------------------------
src          |  production version of PSE code
dev          |  development version of PSE code (may not be stable)
howard       |  test case on Howard's airfoil
parallel     |  test cases on parallel boundary layers
nonpar       |  Nonparallel flow version
notes        |  Code notes and users guide

S. Scott Collis
