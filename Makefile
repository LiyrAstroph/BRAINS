#/*
# * BRAINS
# * (B)LR (R)everberation-mapping Analysis (I)ntegrated with (N)ested (S)ampling
# * Yan-Rong Li, liyanrong@ihep.ac.cn
# * Thu, Aug 4, 2016
# */

SHELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O3 -Wall

#------------target system---------
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
DNEST_INCL  = -I /home/liyropt/Projects/GIT/DNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/DNest -ldnest

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/mbh/mbhd01/soft/gsl/include
GSL_LIBS = -L/mbh/mbhd01/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/mbh/mbhd01/user/liyanrong/soft/mvapich2/lib -lmpich
MPIINCL  = -I/mbh/mbhd01/user/liyanrong/soft/mvapich2/include
LAPACK_INCL = -I/mbh/mbhd01/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/mbh/mbhd01/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
endif

EXEC     = brains
SRC      = ./src
OBJS     = $(SRC)/main.o $(SRC)/allvars.o $(SRC)/read.o $(SRC)/run.o     \
           $(SRC)/dnest_con.o $(SRC)/reconstruct_con.o $(SRC)/init.o     \
           $(SRC)/mathfun.o $(SRC)/dnest_line1d.o $(SRC)/dnest_line2d.o  \
           $(SRC)/reconstruct_line1d.o  $(SRC)/reconstruct_line2d.o      \
           $(SRC)/transfun.o $(SRC)/smooth.o $(SRC)/nrutil.o
 
INCL     = Makefile $(SRC)/allvars.h $(SRC)/proto.h $(SRC)/dnest_con.h   \
           $(SRC)/dnest_line1d.h  $(SRC)/dnest_line2d.h $(SRC)/nrutil.h          

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPICHINCL) $(DNEST_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(MPICHLIB) $(DNEST_LIBS)

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
