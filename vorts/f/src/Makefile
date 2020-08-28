#
#
# note
# * extension .f03 is not supported so I linked the source files to a .f extension
#

# settings for gfortran
FC=gfortran
CFLAGS=-O2 -ffree-form
LFLAGS=

# settings for ifort
# FC=ifort
# CFLAGS=-O3 -fast -free
# LFLAGS=


# create executable
vorts: m_vorts.o main.o
	$(FC) -o ../bin/vorts.exe $(LFLAGS) m_vorts.o main.o

# dependencies
m_vorts.o:./m_vorts.f03
	$(FC) $(CFLAGS) -c ./m_vorts.f03

main.o:./main.f03 m_vorts.o
	$(FC) $(CFLAGS) -c ./main.f03


clean:
	rm -f -r *.o *.mod ../bin/vorts.exe ../bin/vorts
