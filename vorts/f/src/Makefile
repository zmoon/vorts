# This Makefile generates an executable (`../bin/vorts` or `../bin/vorts.exe` on Windows)
# for the main Fortran program.

# Settings for gfortran
FC = gfortran
CFLAGS = -O3 -std=f2008 -Wall -fPIC
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

OBJS2 = \
	m_params.o \
	m_vorts.o \
	capi.o

LIB = ./vorts.so
capi: $(OBJS2)
	$(FC) $(CFLAGS) -shared -o $(LIB) $(LFLAGS) $(OBJS2) 

# Rule to make objects
%.o: %.f90
	$(FC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f -r *.o *.mod $(EXE)
