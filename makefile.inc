CF    = mpif90
FOPT  = -O2 -cpp -dM -fimplicit-none -ffixed-line-length-none -g -Wall -fbacktrace -fcheck=all
FOPT_LINK = -g

FINCL = -I ${FFEMDIR}/mod -I /usr/lib/petsc/include -I /usr/lib/x86_64-linux-gnu/openmpi/include
FLIB = -L${FFEMDIR}/3rdparty/UMFPACK_2.0 -lumfpack_2.0 -lmpi -lmpi_mpifh -lpetsc
