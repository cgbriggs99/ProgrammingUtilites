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

LAPACK_DIR=
