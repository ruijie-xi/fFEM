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

test_moving_exact: $(OBJECT_FILES) test_moving_exact.o $(LIB_FILE)
	${CF} -o $@ $^ $(FLIB) ${FOPT_LINK}

test_static: $(OBJECT_FILES) test_static.o $(LIB_FILE)
	${CF} -o $@ $^ $(FLIB) ${FOPT_LINK}
	
test_mag: $(OBJECT_FILES) test_mag.o $(LIB_FILE)
	${CF} -o $@ $^ $(FLIB) ${FOPT_LINK}
	
clean:
	@rm -f *.o *.mod test_moving_exact test_static test_mag
cleanall:
	@rm -f *.o *.mod output/*.vtk test_moving_exact test_static test_mag
