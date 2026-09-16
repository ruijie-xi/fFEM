.DEFAULT_GOAL := default
.PHONY: all clean lib

all: default

.NOTPARALLEL:
FFEMDIR = ../../..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
	
default: test_multigrid

$(LIBFFEM): lib
$(UMFPACK_LIB): $(LIBFFEM)

module_smoother.o module_prolongation.o module_multigrid.o test_multigrid.o: $(wildcard *.f90) $(LIBFFEM) makefile
module_multigrid.o: module_smoother.o module_prolongation.o
test_multigrid.o: module_multigrid.o

lib: 
	cd ${FFEMDIR} && make && cd -

test_multigrid: module_smoother.o module_prolongation.o module_multigrid.o test_multigrid.o $(LIBFFEM) $(UMFPACK_LIB)
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}
	
clean:
	rm -rf *.o *.mod test_multigrid
