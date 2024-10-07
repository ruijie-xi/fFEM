.PHONY: all clean

FFEMDIR = ..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
%: %.o 
	cd ${FFEMDIR} && make && cd - 
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}
	

test_jacobi_GS: test_jacobi_GS.o module_smoother.o

clean:
	rm -rf *.o *.mod test_jacobi_GS output/ test.txt
