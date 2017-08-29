FILES = constants.f90 routines.f90 switches.f90 variables.f90

all: $(FILES) Makefile
	f2py -m routines -c $(FILES)

clean:
	-rm *.o *.mod

clobber: clean
	-rm -r routines.*.so routines.*.so.dSYM
