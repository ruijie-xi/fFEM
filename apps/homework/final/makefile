.PHONY: all clean lib

FFEMDIR = ../..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
	
default: lib test_multigrid
	
lib: 
	cd ${FFEMDIR} && make && cd -

test_multigrid: module_smoother.o module_prolongation.o module_multigrid.o test_multigrid.o
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}
	
clean:
	rm -rf *.o *.mod test_multigrid
