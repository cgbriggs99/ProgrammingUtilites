include platform.mk

.PHONY:all
all: libconnormath$(LIBEXT) tests docs

libconnormath$(LIBEXT): connormath/libconnormath.so

# Clean up the environment.
.PHONY:clean
clean:
	$(MAKE) -C connormath clean
	$(MAKE) -C docs clean
	$(DEL) *.o *.a *.so *~

# Build the docs.
.PHONY:docs
docs:
	$(MAKE) -C docs html

.PHONY:tests
tests: 
	$(MAKE) -C connormath/tests all

connormath/%:
	$(MAKE) -C connormath $(@F)

