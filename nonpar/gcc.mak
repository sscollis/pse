#==============================================================
#
#  Makefile for pse (GCC)
#
#  Author:  Scott Collis
#
#  Revised: 11-10-97
#
#==============================================================
NAME   = pse
DEBUG  = -g -O2 
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy -c $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -c $(DEBUG)
OFLAGS = $(DEBUG) 
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenBLAS
#LIB    = -lcomplib.sgimath  /usr/people/sashado/lib/bslib.a
FC     = gfortran 
F77    = gfortran
#
#  Define Fortran 90 suffix
#
.SUFFIXES: .f90 
#
#  All modules are listed here
#
MODS = global.o int2str.o fmax.o
#
#  All objects are listed here
#
OBJS = bslib1.o bslib2.o pse.o input.o setup.o solver5.o opper.o genmf.o \
mean.o zeroin.o error.o outf.o plot.o nsolver.o genini.o \
genmean.o genpar.o nonlin.o lst2.o tlst.o adjoint.o tadjoint.o \
asolver_v3.o dasolver_v5.o post.o ddalpha.o

ifdef USE_NR
  ifeq ($(LIBNR_DIR),)
    LIBNR_DIR = $(HOME)/git/NR-utilities
  endif
  LIB += -L$(LIBNR_DIR) -lnr 
endif

#ifdef USE_NR
#  OBJS += nr_rtsafe.o
#endif

$(NAME): $(MODS) $(OBJS) growth.o
	$(FC) $(OFLAGS) $(MODS) $(OBJS) growth.o -o $(NAME) $(LIB)

$(OBJS): $(MODS)

growth.o: fmax.o

blasius: spline.o blasius.o
	$(FC) $(OFLAGS) spline.o blasius.o -o blasius

clean:
	$(RM) *.o *.mod

.f90.o:
	 $(FC) $(F90FLAGS) $*.f90 

.f.o:
	 $(F77) $(FFLAGS) $*.f
