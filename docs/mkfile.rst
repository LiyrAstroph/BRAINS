****************
Makefile
****************

Makefile is used to compile and link the code. There are some important configurations worth mentioning::

  # include spectro-astrometry analysis
  OPTIMIZE += -DSA  
  
  #target system
  SYSTEM="Linux"
  
  ifeq ($(SYSTEM), "Linux")
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

If one wants do spectro-astrometric analysis, switch on the line ``OPTIMIZE += -DSA``.