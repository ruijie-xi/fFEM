

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

FINCL = -I /usr/include/openmpi  \
	-I /usr/lib/petsc/include \
	-I mod

FLIB  = -L/usr/lib -lmpi -lmpi_mpifh \
		-lpetsc

# --- obj file list

OBJECT_FILES = memory_usage.o \
		settings.o \
		mesh_generator.o \
		timer.o \
		writervtk.o \
		solverpetsc.o \
		readmeshfile.o mesh.o\
		quadrature.o \
		fe.o \
		basis.o \
		fe_utils.o \
		visualize.o \
		solver.o \
		assembler.o

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

test_interpolation: ${OBJECTS} $(OBJECT_DIR)/test_interpolation.o
	${CF} -o $(BIN_DIR)/$@ $^ ${FLIB} ${FOPT_LINK}

test_petsc: ${OBJECTS} $(OBJECT_DIR)/test_petsc.o
	${CF} -o $(BIN_DIR)/$@ $^ ${FLIB} ${FOPT_LINK}

test_poisson: ${OBJECTS} $(OBJECT_DIR)/test_poisson.o
	${CF} -o $(BIN_DIR)/$@ $^ ${FLIB} ${FOPT_LINK}

clean:
	@rm $(OBJECT_DIR)/*.o $(MOD_DIR)/*.mod $(BIN_DIR)/*