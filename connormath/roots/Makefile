include ../../platform.mk

CFLAGS+=-I ../include -fPIC

.PHONY: all
all: roots.a

.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~

roots.a: roots.o polyroots.o
	$(AR) -rc $@ $^

roots.o: roots.c

polyroots.o: polyroots.c
polyroots.o: CFLAGS+=--include=$(LAPACK_INCLUDE)
