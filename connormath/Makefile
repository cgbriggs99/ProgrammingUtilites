
include ../platform.mk

%.a:
	$(MAKE) -C $(@D) $(@F)

libconnormath.so: roots/roots.a utils/utils.a
	$(LD) -shared -o $@ --whole-archive $^

.PHONY: clean
clean:
#	$(MAKE) -C funcs clean
#	$(MAKE) -C ints clean
	$(MAKE) -C roots clean
#	$(MAKE) -C stats clean
	$(MAKE) -C tests clean
	$(MAKE) -C utils clean
	$(DEL) *.o *.a *.so *~
