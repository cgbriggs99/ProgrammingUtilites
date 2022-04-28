include ../platform.mk

DEPS=src/roots src/utils
OBJECTS=src/roots/roots.o src/roots/polyroots.o src/utils/utils.o

.PHONY: all
all: libconnormath$(LIBEXT) docs tests

libconnormath$(LIBEXT): $(OBJECTS)
	$(LD) -shared -o $@ $^


# Clean up the environment.
.PHONY: clean
clean:
	$(DEL) *.o *.a *.so *~

# Build the docs.
.PHONY: docs
docs:
	sphinx-build -b html docs/source docs/build/html

.PHONY: tests
tests: src/tests

src/tests: libconnormath$(LIBEXT)
	$(MAKE) -C src/tests all

src/%:
	echo $(@D)
	$(MAKE) -C $(@D) all
