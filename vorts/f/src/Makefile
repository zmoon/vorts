# This Makefile generates an executable (`../bin/vorts.exe`) for the main Fortran program.

# settings for gfortran
FC=gfortran
CFLAGS=-O2 -ffree-form -std=f2003
LFLAGS=

# settings for ifort
# FC=ifort
# CFLAGS=-O3 -fast -free
# LFLAGS=


# create executable
vorts: m_vorts.o main.o
	$(FC) -o ../bin/vorts.exe $(LFLAGS) main.o m_vorts.o m_io.o

# dependencies
m_vorts.o:./m_vorts.f90
	$(FC) $(CFLAGS) -c ./m_vorts.f90

m_io.o:./m_io.f90 m_vorts.o
	$(FC) $(CFLAGS) -c ./m_io.f90

main.o:./main.f90 m_vorts.o m_io.o
	$(FC) $(CFLAGS) -c ./main.f90


clean:
	rm -f -r *.o *.mod ../bin/vorts.exe ../bin/vorts
