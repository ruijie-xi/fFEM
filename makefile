

# --- compiler settings -----------------------------------
SRC_DIR = src
MOD_DIR = mod
OBJECT_DIR = obj
BIN_DIR = bin
VPATH  = $(SRC_DIR):$(MOD_DIR):$(BIN_DIR):$(OBJECT_DIR)

# PETSC_DIR = /home/eureka/petsc
# PETSC_ARCH = arch-linux-c-debug

CF    = gfortran
FOPT  = -O0 -cpp -dM -fimplicit-none -ffixed-line-length-none -g -Wall -fbacktrace -fcheck=all
FOPT_LINK = -g

FINCL = -I mod
FLIB  = -L/usr/lib -L./3rdparty/UMFPACK_2.0 -lumfpack_2.0

# FINCL = $(FINCL) -I /usr/include/openmpi  \
# 	-I /usr/lib/petsc/include \
# FLIB = $(FLIB) -lmpi -lmpi_mpifh -lpetsc

# --- obj file list

OBJECT_FILES = memory_usage.o \
		tools.o \
		quicksort.o \
		settings.o \
		mesh_generator.o \
		timer.o \
		writervtk.o \
		readmeshfile.o \
		mesh.o\
		matvec.o \
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
		solver.o

OBJECTS = $(addprefix $(OBJECT_DIR)/, $(OBJECT_FILES))

# --- commands -------------------------------------------------------

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.F90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.f
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.F
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)


# --- targets --------------------------------------------------------


ffem: ${OBJECTS}
	ar rcs lib/libffem.a $^

test_poisson: ffem $(OBJECT_DIR)/test_poisson.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_magnetic: ffem $(OBJECT_DIR)/test_magnetic.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_matvec: ffem $(OBJECT_DIR)/test_matvec.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

test_quicksort: ffem $(OBJECT_DIR)/test_quicksort.o  
	${CF} -o $(BIN_DIR)/$@ $(word 2,$^) -L./lib -lffem ${FLIB} ${FOPT_LINK}

clean:
	@rm -f $(OBJECT_DIR)/*.o $(MOD_DIR)/*.mod $(BIN_DIR)/* lib/*