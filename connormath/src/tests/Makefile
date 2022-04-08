include ../../../platform.mk

CFLAGS+=-I../../include

test: test.c
test: LDLIBS+= -L ../../ -lconnormath$(LIBEXT) -lm $(LAPACK_LIB)

test_roots.o: test_roots.c
