include ../../../platform.mk

CFLAGS+=-I../../include

all: test.o test_roots run_tests

test.o: test.c

test_roots: test_roots.c test.o
test_roots: CFLAGS+=-L ../../ -lconnormath -lm

run_tests: run_tests.c
