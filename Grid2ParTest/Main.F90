!Entry point to the program

Program main
  use global_mod
  use grid_mod
  use parcel_mod
  use initialise_mod
  use interpolation_mod

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


  print*, "Setting up grids:"

  !now create some grids

  call initialise_grid(ugrid,nxbase, nybase, nzbase)

  call set_grid(ugrid,type=0)

  call initialise_grid(vgrid,nxbase, nybase, nzbase)

  call set_grid(vgrid,type=1)

  call initialise_grid(wgrid,nxbase, nybase, nzbase)

  call set_grid(wgrid,type=2)

  call initialise_grid(pgrid,nxbase*2, nybase*2, nzbase*2)

  call set_grid(pgrid,type=3)

  call initialise_grid(rgrid,nxbase,nybase,nzbase)
  call initialise_grid(sgrid,nxbase,nybase,nzbase)
  call initialise_grid(tgrid,nxbase,nybase,nzbase)



!create some parcels using parcel structures and time their creation
  print *, "Structures:"
  call initialise_parcels(structure=.TRUE.)

  !interpolate grid2par

  call grid2par(grid=ugrid, structure=.true., variable=VEL_X)
  call grid2par(grid=vgrid, structure=.true., variable=VEL_Y)
  call grid2par(grid=wgrid, structure=.true., variable=VEL_Z)
  call par2grid(grid=rgrid, structure=.true.,variable=VORT_X)
  call par2grid(grid=sgrid, structure=.true.,variable=VORT_Y)
  call par2grid(grid=tgrid, structure=.true.,variable=VORT_Z)

  call finalise_parcels(structure=.TRUE.)

  PRINT *, ""

  !create some parcels using arrays and time their creation
  PRINT *, "Arrays:"

  call initialise_parcels(structure=.FALSE.)

  call grid2par(grid=ugrid, structure=.false., variable=VEL_X)
  call grid2par(grid=vgrid, structure=.false., variable=VEL_Y)
  call grid2par(grid=wgrid, structure=.false., variable=VEL_Z)
  call par2grid(grid=rgrid, structure=.false.,variable=VORT_X)
  call par2grid(grid=sgrid, structure=.false.,variable=VORT_Y)
  call par2grid(grid=tgrid, structure=.false.,variable=VORT_Z)

  call finalise_parcels(structure=.FALSE.)


  print*, ""



  !now deallocate grids

  call finalise_grid(ugrid)
  call finalise_grid(vgrid)
  call finalise_grid(wgrid)
  call finalise_grid(pgrid)
  call Finalise_grid(rgrid)
  call Finalise_grid(sgrid)
  call Finalise_grid(tgrid)

end program
