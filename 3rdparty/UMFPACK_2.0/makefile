# Makefile for compiling UMFPACK 2.0

# Compiler and flags
FC = gfortran
FFLAGS = -g -fdefault-real-8 -O2 -std=legacy

# Source and library
SRC_DIR = src
LIB = libumfpack_2.0.a

# List of Fortran source files
SRCS = $(wildcard $(SRC_DIR)/*.f)
OBJS = $(SRCS:.f=.o)

# Default target
all: $(LIB)

# Compile Fortran source files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FFLAGS) -c $< -o $@

# Create the library
$(LIB): $(OBJS)
	ar qc $@ $(OBJS)
	rm -f $(OBJS)

# Clean up
clean:
	rm -f $(LIB)

# Phony targets
.PHONY: all clean