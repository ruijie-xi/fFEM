

# --- compiler settings -----------------------------------
FFEMDIR = .
SRC_DIR = src
MOD_DIR = mod
OBJECT_DIR = obj
BIN_DIR = bin
LIB_DIR = lib
TEST_DIR = test
VPATH  = $(SRC_DIR):$(MOD_DIR):$(BIN_DIR):$(OBJECT_DIR):$(LIB_DIR)

LIBFFEM = libffem.a

include makefile.inc

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


# --- targets --------------------------------------------------------

$(LIBFFEM): ${OBJECTS}
	ar rcs $(LIB_DIR)/$@ $^


clean:
	@rm -f $(OBJECT_DIR)/*.o $(MOD_DIR)/*.mod $(BIN_DIR)/* lib/*.a
	@touch $(BIN_DIR)/.gitkeep
