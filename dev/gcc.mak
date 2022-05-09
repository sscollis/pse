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
#DEBUG  = -O2 -OPT:Olimit=0
DEBUG  = -g -O2
FFLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -std=legacy -c $(DEBUG)
F90FLAGS = -cpp -fdefault-real-8 -fdefault-double-8 -c $(DEBUG)
OFLAGS = $(DEBUG) 
LIB    = -L$(HOME)/local/OpenBLAS/lib -lopenBLAS \
				 -L$(HOME)/local/fftw/lib -lfftw -lrfftw 
FC     = gfortran
F77    = gfortran
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
OBJS = bslib1.o bslib2.o pse.o input.o setup.o opper.o genmf.o mean.o \
zeroin.o error.o outf.o plot.o nsolver.o genini.o \
genmean.o genpar.o nonlin.o lst_v1.o tlst.o adjoint.o tadjoint.o \
asolver_v4.o post.o dasolver_v6.o ddalpha.o duda_v2.o

ifdef USE_NR
	OBJS += nr_rtsafe.o
endif

$(NAME): $(MODS) $(OBJS) growth.o solver_v4.o int2str.o
	$(FC) $(OFLAGS) $(MODS) $(OBJS) solver_v4.o growth.o fmax.o int2str.o \
		-o $(NAME) $(LIB)

$(OBJS): $(MODS)

solver_v4.o: global.o int2str.o fmax.o

asolver_v4.o: global.o int2str.o fmax.o

growth.o: fmax.o

blasius: spline.o blasius.o
	$(FC) $(OFLAGS) spline.o blasius.o -o blasius

clean:
	$(RM) *.o *.mod pse

.f90.o:
	 $(FC) $(F90FLAGS) $*.f90 

.f.o:
	 $(F77) $(FFLAGS) $*.f
