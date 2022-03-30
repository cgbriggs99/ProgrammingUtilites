
# Clean up the environment.
.PHONY:clean
clean:
	python delete_files.py -d *.o *.so *.a *~

# Build the docs.
.PHONY:docs
docs:
	sphinx-build -b html docs/source docs/build/html
