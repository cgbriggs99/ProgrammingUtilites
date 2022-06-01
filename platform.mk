# Contains definitions for platform-dependent options.

ifeq ($(OS),Windows_NT)
LIBEXT=.dll
DEL=del
else
ifeq ($(uname -s),Darwin)
LIBEXT=.lib
else
LDLIBS+=-lm -lpthread
LIBEXT=.so
endif
DEL=rm -rf
endif

LAPACK_INCLUDE=/opt/intel/mkl/include/mkl.h
LAPACK_LIB=-L /opt/intel/mkl/lib/intel64 -lmkl_core
