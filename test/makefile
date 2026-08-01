.PHONY: all clean

FFEMDIR = ..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
%: %.o 
	cd ${FFEMDIR} && make && cd - 
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}
	
executables = test_mesh test_interpolation test_poisson test_stokes test_matvec test_quicksort test_mpi test_poisson_parallel
	
all: ${executables}

test_mesh: test_mesh.o 

test_interpolation: test_interpolation.o 

test_poisson: test_poisson.o 

test_stokes: test_stokes.o 

test_matvec: test_matvec.o 

test_quicksort: test_quicksort.o 

test_mpi: test_mpi.o 

clean: 
	rm -rf *.o *.mod ${executables} output/ test.txt
