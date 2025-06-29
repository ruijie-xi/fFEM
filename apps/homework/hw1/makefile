.PHONY: all clean lib

FFEMDIR = ../..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
	
default: lib test_jacobi_GS
	
lib: 
	cd ${FFEMDIR} && make && cd -

test_jacobi_GS: module_smoother.o test_jacobi_GS.o 
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}

clean:
	rm -rf *.o *.mod test_jacobi_GS output/ test.txt *.png
