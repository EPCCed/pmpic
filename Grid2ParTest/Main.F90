!Entry point to the program

Program main
  use global_mod
  use grid_mod
  use parcel_mod
  use initialise_mod
  use interpolation_mod

  implicit none

  integer :: i
  integer, parameter :: numits=20


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

!$OMP parallel
!$OMP MASTER
  print *, "OMP_NUM_THREADS=",OMP_GET_NUM_THREADS()
  print *, ""
!$OMP END MASTER
!$OMP END PARALLEL


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

print*, ""

!create some parcels using parcel structures and time their creation
  print *, "Structures:"
  call initialise_parcels(structure=.TRUE.,shuffle=.True.)

  !interpolate grid2par

  tp2g=0.d0
  tg2p=0.d0

  do i=1,numits
      call grid2par(grid=ugrid, structure=.true., variable=VEL_X)
      call grid2par(grid=vgrid, structure=.true., variable=VEL_Y)
      call grid2par(grid=wgrid, structure=.true., variable=VEL_Z)
      call par2grid(grid=rgrid, structure=.true.,variable=VORT_X)
      call par2grid(grid=sgrid, structure=.true.,variable=VORT_Y)
      call par2grid(grid=tgrid, structure=.true.,variable=VORT_Z)
  enddo

  print *, "mean grid2par time=",tg2p/numits/3.
  print *, "mean par2grid time=",tp2g/numits/3.


  call finalise_parcels(structure=.TRUE.)

  PRINT *, ""

  !create some parcels using arrays and time their creation
  PRINT *, "Arrays:"

  call initialise_parcels(structure=.FALSE., shuffle=.true.)

  tp2g=0.d0
  tg2p=0.d0

  do i=1,numits
      call grid2par(grid=ugrid, structure=.false., variable=VEL_X)
      call grid2par(grid=vgrid, structure=.false., variable=VEL_Y)
      call grid2par(grid=wgrid, structure=.false., variable=VEL_Z)
      call par2grid(grid=rgrid, structure=.false.,variable=VORT_X)
      call par2grid(grid=sgrid, structure=.false.,variable=VORT_Y)
      call par2grid(grid=tgrid, structure=.false.,variable=VORT_Z)
  enddo

  print *, "mean grid2par time=",numits/15./3.
  print *, "mean par2grid time=",numits/15./3.

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
