include ../platform.mk

.PHONY:all
all: library docs

.PHONY:library
library: src/all

# Clean up the environment.
.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~

# Build the docs.
.PHONY:docs
docs:
	sphinx-build -b html docs/source docs/build/html
