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

  call finalise_parcels(structure=.TRUE.)

  PRINT *, ""

  !create some parcels using arrays and time their creation
  PRINT *, "Arrays:"

  call initialise_parcels(structure=.FALSE.)

  call finalise_parcels(structure=.FALSE.)


  !now create some grids

  call initialise_grid(ugrid,nxbase, nybase, nzbase)

  call set_grid(ugrid,type=0)

  call initialise_grid(vgrid,nxbase, nybase, nzbase)

  call set_grid(vgrid,type=1)

  call initialise_grid(wgrid,nxbase, nybase, nzbase)

  call set_grid(wgrid,type=2)

  call initialise_grid(pgrid,nxbase, nybase, nzbase)

  call set_grid(pgrid,type=3)

  !now deallocate grids

  call finalise_grid(ugrid)
  call finalise_grid(vgrid)
  call finalise_grid(wgrid)
  call finalise_grid(pgrid)

end program
