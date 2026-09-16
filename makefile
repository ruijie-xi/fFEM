

.DEFAULT_GOAL := all
.NOTPARALLEL:
.PHONY: all default clean clean-all clean-umfpack check

# --- compiler settings -----------------------------------
FFEMDIR = .
SRC_DIR = src
MOD_DIR = mod
OBJECT_DIR = obj
BIN_DIR = bin
LIB_DIR = lib
TEST_DIR = test
VPATH  = $(SRC_DIR):$(MOD_DIR):$(BIN_DIR):$(OBJECT_DIR):$(LIB_DIR)

LIBFFEM = $(LIB_DIR)/libffem.a

include makefile.inc

# --- obj file list ----------------------------------------
OBJECT_FILES = memory_usage.o \
		timer.o \
		tools.o \
		writervtk.o \
		readmeshfile.o \
		quicksort.o \
		settings.o \
		geometry.o \
		matvec.o \
		mesh_generator.o \
		mesh.o\
		quadrature.o \
		fespace_P0.o \
		fespace_DG1.o \
		fespace_P1.o \
		fespace_Q0.o \
		fespace_Q1.o \
		fespace_Q2.o \
		fespace_P2.o \
		fespace_QuadNedelec1.o \
		fespace_QuadRT2.o \
		fespace_QuadRT1.o \
		fe.o \
		fe_utils.o \
		visualize.o \
		assembler.o \
		solver_umfpack2.o
		
ifeq ($(USE_PETSC),1)
OBJECT_FILES += solver_petsc.o
endif

OBJECTS = $(addprefix $(OBJECT_DIR)/, $(OBJECT_FILES))

# Conservative module dependency handling for this small, ordered build.
$(OBJECTS): $(wildcard $(SRC_DIR)/*.f90) makefile makefile.inc

# --- commands -------------------------------------------------------

$(OBJECT_DIR)/%.o: $(SRC_DIR)/%.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ -J $(MOD_DIR)


# --- targets --------------------------------------------------------

default: $(LIBFFEM)

all: $(UMFPACK_LIB) $(LIBFFEM)

$(UMFPACK_LIB): $(wildcard $(UMFPACK_DIR)/src/*.f) $(UMFPACK_DIR)/makefile
	$(MAKE) -C $(UMFPACK_DIR)

$(LIBFFEM): ${OBJECTS}
	ar rcs $@ $^

check: all
	$(MAKE) -C test check

clean-umfpack:
	cd $(UMFPACK_DIR) && make clean && cd -

clean:
	@rm -f $(OBJECT_DIR)/*.o $(MOD_DIR)/*.mod $(BIN_DIR)/* lib/*.a
	@touch $(BIN_DIR)/.gitkeep

clean-all: clean clean-umfpack
