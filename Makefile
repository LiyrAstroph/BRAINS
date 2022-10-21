#/*
# * BRAINS
# * (B)LR (R)everberation-mapping (A)nalysis (I)n AGNs with (N)ested (S)ampling
# * Yan-Rong Li, liyanrong@ihep.ac.cn
# * Thu, Aug 4, 2016
# */

SHELL=/bin/bash

CC       = mpicc
OPTIMIZE = -O2 -Wall -finline-functions -fcommon
#OPTIMIZE += -DDebug

# include spectro-astrometry analysis
#OPTIMIZE += -DSpecAstro

# get GIT description
GITCHECK := $(shell git 2>/dev/null)
ifneq ($(strip $(GITCHECK)),)
	GITREPO := $(shell git rev-parse --is-inside-work-tree 2>/dev/null)
	ifneq ($(findstring true, $(GITREPO)),)
		GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
		OPTIMIZE += -DGITVERSION=\"$(GIT_VERSION)\"
	    
		GIT_DATE := $(firstword $(shell git --no-pager show --date=short --format="%ad" --name-only))
		OPTIMIZE += -DGITDATE=\"$(GIT_DATE)\"
	endif
endif

#------------target system---------
#SYSTEM="Darwin"
SYSTEM="Linux"
#SYSTEM="Cluster"
#SYSTEM="TianheII"

ifeq ($(SYSTEM), "Linux")
NCORE      :=$(grep -c ^processor /proc/cpuinfo)
GSL_INCL    = $(shell pkg-config --cflags gsl) 
GSL_LIBS    = $(shell pkg-config --libs gsl) 
LAPACK_INCL = -I/usr/include/lapacke
LAPACK_LIBS = -L/usr/lib64 -llapacke -llapack -lblas
DNEST_INCL  = -I /home/liyropt/Projects/GIT/CDNest/
DNEST_LIBS  = -L /home/liyropt/Projects/GIT/CDNest -ldnest
FFTW_INCL   = $(shell pkg-config --cflags fftw3) 
FFTW_LIBS   = $(shell pkg-config --libs fftw3) 

MPICHINCL     = $(shell pkg-config --cflags mpich) 
MPICHLIB    = $(shell pkg-config --libs mpich) 
endif

ifeq ($(SYSTEM), "Cluster")
GSL_INCL = -I/sharefs/mbh/user/liyanrong/soft/gsl/include
GSL_LIBS = -L/sharefs/mbh/user/liyanrong/soft/gsl/lib  -lgsl -lgslcblas -lm
MPICHLIB = -L/sharefs/mbh/user/liyanrong/soft/mpich3/lib -lmpich
MPIINCL  = -I/sharefs/mbh/user/liyanrong/soft/mpich3/include
LAPACK_INCL = -I/sharefs/mbh/user/liyanrong/soft/lapack/include
LAPACK_LIBS = -L/sharefs/mbh/user/liyanrong/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
FFTW_INCL = -I/sharefs/mbh/user/liyanrong/soft/fftw/include
FFTW_LIBS = -L/sharefs/mbh/user/liyanrong/soft/fftw/lib -lfftw3

DNEST_INCL  = -I /sharefs/mbh/user/liyanrong/GIT/DNest/
DNEST_LIBS  = -L /sharefs/mbh/user/liyanrong/GIT/DNest -ldnest
endif

ifeq ($(SYSTEM), "TianheII")
GSL_INCL =
GSL_LIBS = -lgsl -lgslcblas -lm
MPICHLIB = -lmpich
MPIINCL  =
LAPACK_INCL = -I/HOME/ihep_yrli_1/BIGDATA/soft/lapack/include
LAPACK_LIBS = -L/HOME/ihep_yrli_1/BIGDATA/soft/lapack/lib -llapacke -llapack -lblas -lgfortran
FFTW_INCL =
FFTW_LIBS = -lfftw3

DNEST_INCL  = -I /HOME/ihep_yrli_1/BIGDATA/soft/DNest/
DNEST_LIBS  = -L /HOME/ihep_yrli_1/BIGDATA/soft/DNest -ldnest
endif


EXEC     = brains
SRC      = ./src
OBJS     = $(SRC)/main.o $(SRC)/allvars.o $(SRC)/read.o $(SRC)/run.o     \
           $(SRC)/dnest_con.o $(SRC)/reconstruct_con.o $(SRC)/init.o     \
           $(SRC)/mathfun.o $(SRC)/dnest_line1d.o $(SRC)/dnest_line2d.o  \
           $(SRC)/reconstruct_line1d.o  $(SRC)/reconstruct_line2d.o      \
           $(SRC)/transfun.o $(SRC)/smooth_fftw.o $(SRC)/nrutil.o        \
           $(SRC)/system.o  $(SRC)/sim.o $(SRC)/help.o $(SRC)/version.o  \
           $(SRC)/blr_models.o $(SRC)/command_line.o                     \
           $(SRC)/user_blr_model.o   $(SRC)/user_transfun.o              \
           $(SRC)/reconstruct_sa.o   $(SRC)/dnest_sa.o                   \
           $(SRC)/specastro.o $(SRC)/dnest_sa1d.o $(SRC)/dnest_sa2d.o    \
           $(SRC)/reconstruct_sa1d.o  $(SRC)/reconstruct_sa2d.o          \
           $(SRC)/sarm.o  $(SRC)/reconstruct_sarm.o  $(SRC)/dnest_sarm.o \
           $(SRC)/sa_gravity.o          
 
INCL     = Makefile $(SRC)/allvars.h $(SRC)/proto.h $(SRC)/dnest_con.h   \
           $(SRC)/dnest_line1d.h  $(SRC)/dnest_line2d.h $(SRC)/nrutil.h  \
           $(SRC)/blr_models.h   $(SRC)/brains.h $(SRC)/version.h        \
           $(SRC)/command_line.h  $(SRC)/user_blr_model.h                \
           $(SRC)/user_transfun.h  $(SRC)/dnest_sa.h $(SRC)/dnest_sa1d.h \
           $(SRC)/dnest_sa2d.h  $(SRC)/mathfun.h  $(SRC)/dnest_sarm.h    \
           $(SRC)/sa_gravity.h

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) $(GSL_INCL) $(LAPACK_INCL) $(MPICHINCL) $(DNEST_INCL) $(FFTW_INCL)
LIBS     = $(GSL_LIBS) $(LAPACK_LIBS) $(MPICHLIB) $(DNEST_LIBS) $(FFTW_LIBS)

$(EXEC):$(OBJS)
	cd $(SRC)
	$(CC) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS) $(EXEC)
