.DEFAULT_GOAL := default
.PHONY: all clean lib

all: default

.NOTPARALLEL:
FFEMDIR = ../../..
LIBFFEM = ${FFEMDIR}/lib/libffem.a

include ${FFEMDIR}/makefile.inc

%.o: %.f90
	${CF} -c ${FOPT} ${FINCL} $< -o $@ 
	
	
default: test_jacobi_GS

$(LIBFFEM): lib
$(UMFPACK_LIB): $(LIBFFEM)

module_smoother.o test_jacobi_GS.o: $(wildcard *.f90) $(LIBFFEM) makefile
test_jacobi_GS.o: module_smoother.o

lib: 
	cd ${FFEMDIR} && make && cd -

test_jacobi_GS: module_smoother.o test_jacobi_GS.o $(LIBFFEM) $(UMFPACK_LIB)
	${CF} -o $@ $^ -L${FFEMDIR}/lib -lffem ${FLIB} ${FOPT_LINK}

clean:
	rm -rf *.o *.mod test_jacobi_GS output/ test.txt *.png
