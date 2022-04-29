include ../../platform.mk

all: utils.a

utils.a: utils.o
	$(AR) -rc $@ $^

utils.o: utils.c
utils.o: CFLAGS+=-I../include

.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~
