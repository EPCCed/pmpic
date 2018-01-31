FC=mpif90
FFLAGS= -O2 -Wall -fbounds-check -fopenmp

FFILES=Main.f90 grid.o parcel.o initialise.o global.o

%.mod: %.F90
	$(FC) $(FFLAGS) -c $<

main : Main.f90 global.mod initialise.mod global.mod parcel.mod grid.mod
	$(FC) $(FFLAGS) $(FFILES) -o main

global.mod : global.F90

grid.mod : grid.F90 global.mod

parcel.mod : parcel.F90 global.mod

initialise.mod : initialise.F90 global.mod parcel.mod grid.mod


clean:
	rm *.o *.mod
