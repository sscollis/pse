#==============================================================
#
#  Makefile for pse (Darwin)
#
#  Author:  Scott Collis
#
#  Revised: 12-23-2019 
#
#==============================================================
NAME   = pse
DEBUG  = -g -O2
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy -c $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -c $(DEBUG)
OFLAGS = -fdefault-real-8 -fdefault-double-8 $(DEBUG) 
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenBLAS
#LIB   = -lcomplib.sgimath /usr/people/collis/lib/bslib.a
FC     = gfortran
F77    = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = global.o int2str.o
#
#  All objects are listed here
#
OBJS = bslib1.o bslib2.o pse.o input.o setup.o opper.o genmf.o mean.o \
zeroin.o error.o outf.o plot.o nsolver.o genini.o \
genmean.o genpar.o nonlin.o lst_v1.o tlst.o adjoint.o tadjoint.o \
asolver_v4.o post.o dasolver_v6.o ddalpha.o duda_v2.o

OBJS2 = growth.o solver_v4.o fmax.o

ifdef USE_NR
	OBJS += nr_rtsafe.o
endif

$(NAME): $(MODS) $(OBJS) $(OBJS2)
	$(FC) $(OFLAGS) $(MODS) $(OBJS) $(OBJS2) $(LIB) -o $(NAME)

$(OBJS): $(MODS)

growth.o: fmax.o

solver_v4.o: global.o fmax.o

blasius: spline.o blasius.o
	$(FC) $(OFLAGS) spline.o blasius.o -o blasius

clean:
	$(RM) *.o *.mod $(NAME)

.f90.o:
	 $(FC) $(F90FLAGS) $*.f90 

.f.o:
	 $(F77) $(FFLAGS) $*.f
