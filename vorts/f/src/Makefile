# This Makefile generates an executable (`../bin/vorts.exe`) for the main Fortran program.

# settings for gfortran
FC = gfortran
CFLAGS = -O3 -std=f2008
LFLAGS =

# settings for ifort
# FC = ifort
# CFLAGS = -O3 -fast -free
# LFLAGS =


# Create executable
OBJS = main.o m_vorts.o m_io.o
vorts: $(OBJS)
	$(FC) $(CFLAGS) -o ../bin/vorts.exe $(LFLAGS) $(OBJS)

# Objects
m_vorts.o:./m_vorts.f90
	$(FC) $(CFLAGS) -c $<

m_io.o:./m_io.f90 m_vorts.o
	$(FC) $(CFLAGS) -c $<

main.o:./main.f90 m_vorts.o m_io.o
	$(FC) $(CFLAGS) -c $<

.PHONY: clean
clean:
	rm -f -r *.o *.mod ../bin/vorts.exe ../bin/vorts
