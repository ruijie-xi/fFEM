CF    = mpif90
FOPT  = -O2 -cpp -dM -fimplicit-none -ffixed-line-length-none -g -Wall -fbacktrace -fcheck=all
FOPT_LINK = -g

USE_PETSC = 0
PETSC_DIR = /usr/lib/petsc

FINCL = -I ${FFEMDIR}/mod
FLIB = -L${FFEMDIR}/3rdparty/UMFPACK_2.0 -lumfpack_2.0 

ifeq ($(USE_PETSC),1)
  FINCL += -I${PETSC_DIR}/include
  FLIB += -L${PETSC_DIR}/lib -lpetsc
endif
