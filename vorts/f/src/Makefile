# This Makefile generates an executable (`../bin/vorts` or `../bin/vorts.exe` on Windows)
# for the main Fortran program.

# Settings for gfortran
FC = gfortran
CFLAGS = -O3 -std=f2008
LFLAGS =

# Settings for ifort
# FC = ifort
# CFLAGS = -O3 -fast -free
# LFLAGS =

# Objects (order matters)
OBJS = \
	m_params.o \
	m_vorts.o \
	m_io.o \
	main.o

# Create executable
EXE = ../bin/vorts
vorts: $(OBJS)
	$(FC) $(CFLAGS) -o $(EXE) $(LFLAGS) $(OBJS)

# Rule to make objects
.SUFFIXES: .f90
.f90.o:
	$(FC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f -r *.o *.mod $(EXE)
