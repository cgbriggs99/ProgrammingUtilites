include ../platform.mk

DEPS=src/roots
OBJECTS=src/roots/roots.o src/roots/polyroots.o

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
tests: src/tests/test_roots

src/tests/%:
	$(MAKE) -C src/tests $(@F)

src/%:
	$(MAKE) -C src $(@F)
