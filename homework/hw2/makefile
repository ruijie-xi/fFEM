.PHONY: all clean lib

FFEMDIR = ../..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
	
default: lib test_twogrid
	
lib: 
	cd ${FFEMDIR} && make && cd -

test_twogrid: module_smoother.o test_twogrid.o module_prolongation.o
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}

clean:
	rm -rf *.o *.mod test_twogrid output/ test.txt *.png
