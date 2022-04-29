include ../../platform.mk

CFLAGS+=-I../include

all: test_roots test_utils

test_roots: test_roots.o
test_roots: LDLIBS+=-L ../ -lconnormath -lm $(LAPACK_LIB)
test_roots: LDFLAGS+=-Wl,-R$(abspath ../)

test_utils: test_utils.o
test_utils: LDLIBS+=-L ../ -lconnormath -lm $(LAPACK_LIB)
test_utils: LDFLAGS+=-Wl,-R$(abspath ../)

test_roots.o: test_roots.c

test_utils.o: test_utils.c

.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~ test_roots
