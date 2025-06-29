FFEM_DIR = ../..
# --- compiler settings -----------------------------------

CF    = gfortran
FOPT  = -O -fimplicit-none -g -Wall -fbacktrace -fcheck=all -std=gnu
FOPT_LINK = -g

FINCL = -I$(FFEM_DIR)/mod 
FLIB = -L$(FFEM_DIR)/lib -lffem -L$(FFEM_DIR)/3rdparty/UMFPACK_2.0 -lumfpack_2.0

# --- obj file list ----------------------------------------
OBJECT_FILES = hydrodynamic.o magnetic.o

LIB_FILE = $(FFEM_DIR)/lib/libffem.a

%.o: %.f90
	$(CF) -c $(FOPT) $(FINCL) $< -o $@

%.o: %.f
	$(CF) -c $(FOPT) $(FINCL) $< -o $@ 


# --- targets --------------------------------------------------------

magnetic_module: magnetic.o
	
clean:
	@rm -f *.o *.mod 
cleanall:
	@rm -f *.o *.mod output/*.vtk
