
ifeq ($(OS),Windows_NT)
DEL=del
else
DEL=rm -rf
endif

# Clean up the environment.
.PHONY:clean
clean:
	$(DEL) *.o *.a *.so *~

# Build the docs.
.PHONY:docs
docs:
	sphinx-build -b html docs/source docs/build/html
