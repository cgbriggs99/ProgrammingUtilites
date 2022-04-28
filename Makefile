
ifeq ($(OS),Windows_NT)
LIBEXT=.dll
DEL=del
else
ifeq ($(uname -s),Darwin)
LIBEXT=.lib
else
LIBEXT=.so
endif
DEL=rm -rf
endif

all: 

# Clean up the environment.
.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~

# Build the docs.
.PHONY:docs
docs:
	sphinx-build -b html docs/source docs/build/html
