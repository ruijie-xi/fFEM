CF    = mpif90
FOPT  = -O -cpp -dM -fimplicit-none -ffixed-line-length-none -g -Wall -fbacktrace -fcheck=all
FOPT_LINK = -g

UMFPACK_DIR = $(FFEMDIR)/3rdparty/UMFPACK_2.0
UMFPACK_LIB = $(UMFPACK_DIR)/libumfpack_2.0.a
USE_PETSC = 0
PETSC_DIR = /usr/lib/petsc

FINCL = -I ${FFEMDIR}/mod
FLIB = -L${UMFPACK_DIR} -lumfpack_2.0 

ifeq ($(USE_PETSC),1)
  FINCL += -I${PETSC_DIR}/include
  FLIB += -L${PETSC_DIR}/lib -lpetsc
endif
