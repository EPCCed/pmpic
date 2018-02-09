FC=mpif90
FFLAGS= -O2 -Wall -fbounds-check -fopenmp

FFILES=main.o grid.o parcel.o initialise.o global.o interpolation.o

%.o: %.F90
	$(FC) $(FFLAGS) -c $<

main : $(FFILES)
	$(FC) $(FFLAGS) $(FFILES) -o main

main.o: main.F90 global.o grid.o parcel.o initialise.o interpolation.o

grid.o : grid.F90 global.o

parcel.o : parcel.F90 global.o

initialise.o : initialise.F90 global.o parcel.o grid.o

interpolation.o : global.o parcel.o grid.o


clean:
	rm *.o *.mod
