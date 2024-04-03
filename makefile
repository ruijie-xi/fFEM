

# --- compiler settings -----------------------------------
SRC_DIR = src
MOD_DIR = mod
OBJECT_DIR = obj
BIN_DIR = bin
LIB_DIR = lib
TEST_DIR = test
VPATH  = $(SRC_DIR):$(MOD_DIR):$(BIN_DIR):$(OBJECT_DIR):$(LIB_DIR)

LIBFFEM = libffem.a

# PETSC_DIR = /home/eureka/petsc
# PETSC_ARCH = arch-linux-c-debug

CF    = mpif90
FOPT  = -O2 -cpp -dM -fimplicit-none -ffixed-line-length-none -g -Wall -fbacktrace -fcheck=all
FOPT_LINK = -g

FINCL = -I mod -I /usr/lib/petsc/include -I /usr/lib/x86_64-linux-gnu/openmpi/include
FLIB = -L./3rdparty/UMFPACK_2.0 -lumfpack_2.0 -lmpi -lmpi_mpifh -lpetsc

# --- obj file list ----------------------------------------
OBJECT_FILES = memory_usage.o \
		timer.o \
		tools.o \
		writervtk.o \
		readmeshfile.o \
		quicksort.o \
		settings.o \
		matvec.o \
		mesh_generator.o \
		mesh.o\
		quadrature.o \
		fespace_P0.o \
		fespace_DG1.o \
		fespace_P1.o \
		fespace_Q0.o \
		fespace_Q1.o \
		fespace_P2.o \
		fespace_QuadNedelec1.o \
		fespace_QuadRT1.o \
		fe.o \
		fe_utils.o \
		visualize.o \
		assembler.o \
		solver_umfpack2.o \
		solver_petsc.o

OBJECTS = $(addprefix $(OBJECT_DIR)/, $(OBJECT_FILES))

# --- commands -------------------------------------------------------

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)

$(OBJECT_DIR)/%.o: $(TEST_DIR)/%.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)


# --- targets --------------------------------------------------------

$(LIBFFEM): ${OBJECTS}
	ar rcs $(LIB_DIR)/$@ $^

test_mesh: $(LIBFFEM) $(OBJECT_DIR)/test_mesh.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_interpolation: $(LIBFFEM) $(OBJECT_DIR)/test_interpolation.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_poisson: $(LIBFFEM) $(OBJECT_DIR)/test_poisson.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_stokes: $(LIBFFEM) $(OBJECT_DIR)/test_stokes.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_magnetic: $(LIBFFEM) $(OBJECT_DIR)/test_magnetic.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_matvec: $(LIBFFEM) $(OBJECT_DIR)/test_matvec.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_quicksort: $(LIBFFEM) $(OBJECT_DIR)/test_quicksort.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_mpi: $(LIBFFEM) $(OBJECT_DIR)/test_mpi.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_poisson_parallel: $(LIBFFEM) $(OBJECT_DIR)/test_poisson_parallel.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

clean:
	@rm -f $(OBJECT_DIR)/*.o $(MOD_DIR)/*.mod $(BIN_DIR)/* lib/*.a
	@touch $(BIN_DIR)/.gitkeep
