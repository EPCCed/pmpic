FC=mpif90
FFLAGS= -O2 -Wall -fbounds-check -fopenmp

FC=ftn
FFLAGS=-h omp
#FFLAGS=-O2 -openmp
#FFLAGS=-O2 -fopenmp

FFILES=Main.o grid.o parcel.o initialise.o global.o interpolation.o testInterp.o

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

main : $(FFILES)
	$(FC) $(FFLAGS) $(FFILES) -o main

Main.o: Main.F90 global.o grid.o parcel.o initialise.o interpolation.o testInterp.o

grid.o : grid.F90 global.o

parcel.o : parcel.F90 global.o

initialise.o : initialise.F90 global.o parcel.o grid.o

interpolation.o : global.o parcel.o grid.o

testInterp.o : global.o parcel.o grid.o interpolation.o


clean:
	rm *.o *.mod
