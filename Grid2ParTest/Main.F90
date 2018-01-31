!Entry point to the program

Program main
  use global_mod
  use grid_mod
  use parcel_mod
  use initialise_mod

  implicit none


  !setup computational domain

  xmin=-3.d0
  xmax=3.d0
  ymin=-3.d0
  ymax=3.d0
  zmin=0.d0
  zmax=3.d0

  nxbase=128
  nybase=128
  nzbase=64

!create some parcels using parcel structures and time their creation
  print *, "Structures:"
  call initialise_parcels(structure=.TRUE.)

  call finalize_parcels(structure=.TRUE.)

  PRINT *, ""

  !create some parcels using arrays and time their creation
  PRINT *, "Arrays:"

  call initialise_parcels(structure=.FALSE.)

  call finalize_parcels(structure=.FALSE.)

end program
