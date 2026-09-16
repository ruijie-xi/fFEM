.DEFAULT_GOAL := all
.NOTPARALLEL:
.PHONY: all check clean FORCE

FFEMDIR = ..
LIBFFEM = $(FFEMDIR)/lib/libffem.a
include $(FFEMDIR)/makefile.inc

executables = test_mesh test_interpolation test_poisson test_stokes test_matvec test_quicksort test_mpi test_regression

all: $(executables)

$(LIBFFEM): FORCE
	$(MAKE) -C $(FFEMDIR) all

$(executables): %: %.o $(LIBFFEM) $(UMFPACK_LIB)
	$(CF) -o $@ $< -L$(FFEMDIR)/lib -lffem $(FLIB) $(FOPT_LINK)

$(UMFPACK_LIB): $(LIBFFEM)

%.o: %.f90 $(LIBFFEM) $(FFEMDIR)/makefile.inc
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

check: all
	python3 run_tests.py

clean:
	rm -f *.o *.mod $(executables)
