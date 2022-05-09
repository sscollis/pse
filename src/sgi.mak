#==============================================================
#
#  Makefile for pse (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 11-10-97
#
#==============================================================
NAME   = pse
DEBUG  = -O2 -OPT:Olimit=0
#DEBUG  = -g
FFLAGS = -cpp -r8 -c $(DEBUG)
F90FLAGS = -cpp -r8 -c $(DEBUG)
OFLAGS = -r8 $(DEBUG) 
LIB    = -lcomplib.sgimath /usr/people/collis/lib/bslib.a
COMP   = f90
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = global.o 
#
#  All objects are listed here
#
OBJS = \
pse.o input.o setup.o opper.o genmf.o mean.o \
zeroin.o error.o outf.o plot.o nsolver.o genini.o \
genmean.o genpar.o nonlin.o lst_v1.o tlst.o adjoint.o tadjoint.o \
asolver_v4.o rtsafe.o post.o dasolver_v6.o ddalpha.o duda_v2.o

OBJS2 = growth.o solver_v4.o int2str.o fmax.o

$(NAME): $(MODS) $(OBJS) $(OBJS2)
	$(COMP) $(OFLAGS) $(MODS) $(OBJS) $(OBJS2) -o $(NAME) $(LIB)

$(OBJS): $(MODS)

growth.o: fmax.o

solver_v4.o: global.o int2str.o fmax.o

blasius: spline.o blasius.o
	$(COMP) $(OFLAGS) spline.o blasius.o -o blasius

clean:
	/bin/rm *.o *.mod

.f90.o:
	 $(COMP) $(F90FLAGS) $*.f90 

.f.o:
	 f77 -col120 $(FFLAGS) $*.f
