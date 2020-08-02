#
#
# this is an ifort version
# extension .f03 is not supported so I linked the source files to a .f extension
#

# settings for gfortran
FC=gfortran 
CFLAGS=-O2
LFLAGS=

# settings for ifort
# FC=ifort
# CFLAGS=-O3 -fast -free
# LFLAGS=


# create executable
vorts: m_vorts.o main.o
	$(FC) -o ../bin/vorts.exe $(LFLAGS) m_vorts.o main.o

# dependencies
m_vorts.o:./m_vorts.f
	$(FC) $(CFLAGS) -c ./m_vorts.f

main.o:./main.f m_vorts.o
	$(FC) $(CFLAGS) -c ./main.f


clean:
	rm -f -r *.o *.mod ../bin/vorts.exe ../bin/vorts


