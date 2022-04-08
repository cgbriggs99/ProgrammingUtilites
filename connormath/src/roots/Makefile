include ../../../platform.mk

CFLAGS+=-I ../../include -fPIC

.PHONY: all
all: roots.o polyroots.o

roots.o: roots.c

polyroots.o: polyroots.c
polyroots.o: CFLAGS+=--include=$(LAPACK_INCLUDE)
