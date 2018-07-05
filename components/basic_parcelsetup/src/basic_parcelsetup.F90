!basic parcel setup - parcels are created to fill the volume with n_per_dir parcels per
!grid cell per direction. Only x, y and z and vol. All the other proeprties are zero

module basic_parcelsetup_mod
  use datadefn_mod, only : DEFAULT_PRECISION, PARCEL_INTEGER, MPI_PARCEL_INT
  use state_mod, only: model_state_type
  use optionsdatabase_mod, only: options_get_integer
  use parcel_interpolation_mod, only: x_coords, y_coords, z_coords
  use MPI
  use timer_mod
  use monc_component_mod, only: component_descriptor_type

contains

  type(component_descriptor_type) function basic_parcelsetup_get_descriptor()
    basic_parcelsetup_get_descriptor%name="basic_parcelsetup"
    basic_parcelsetup_get_descriptor%version=0.1
    basic_parcelsetup_get_descriptor%initialisation=>initialisation_callback

  end function basic_parcelsetup_get_descriptor

  subroutine initialisation_callback(state)
    type(model_state_type), intent(inout), target :: state
    integer :: n_per_dir
    integer :: n_per_cell
    integer :: nx, ny, nz
    real(kind=DEFAULT_PRECISION) :: dx, dy, dz
    integer :: nnx, nny, nnz
    integer(kind=PARCEL_INTEGER) :: nparcels
    integer(kind=PARCEL_INTEGER) :: n
    integer :: i, j, k
    integer :: xstart, xstop, ystart, ystop, zstart, zstop
    real(kind=DEFAULT_PRECISION) :: ddx, ddy, ddz, xs, ys, zs
    integer :: ii, jj, kk
    integer :: q
    integer :: handle

    call register_routine_for_timing("Basic_Setup",handle,state)
    call timer_start(handle)

    nx = state%local_grid%size(3) + 2*state%local_grid%halo_size(3)
    ny = state%local_grid%size(2) + 2*state%local_grid%halo_size(2)
    nz = state%local_grid%size(1) + 2*state%local_grid%halo_size(1)

    n_per_dir=options_get_integer(state%options_database,"parcels_per_cell_dir")



    !calculate number of parcels in a cell
    n_per_cell=n_per_dir**3

    !calculate number of cells in each direction (total grid size minus halo cells)
    nnx=nx - 2*state%local_grid%halo_size(3) !- 1
    nny=ny - 2*state%local_grid%halo_size(2) !- 1
    nnz=nz - 2*state%local_grid%halo_size(1) - 1

    nparcels=n_per_cell*(nnx)*(nny)*(nnz)

    if (state%parallel%my_rank .eq. 0) write(*,"(i2.1,a)")  n_per_dir, " parcels per direction per cell"
    if (state%parallel%my_rank .eq. 0) write(*,"(a,i4.1,a)") " So ", n_per_cell," Parcels per cell"
    !if (state%parallel%my_rank .eq. 0) write(*,"(a,i12,a,i12)") "Setting up", nparcels,"parcels in", nnx*nny*nnz, "cells"

    !print*, "nnx, nny, nnz=", nnx, nny, nnz
    !print*, "maxparcels=",state%parcels%maxparcels_local,  " nparcels=", nparcels

    if (nparcels .gt. state%parcels%maxparcels_local) then
      error stop "Maxparcels is not big enough for the number of parcels per cell requested"
    endif

    state%parcels%numparcels_local = nparcels

    !set global parcel count
    call MPI_Allreduce(sendbuf=state%parcels%numparcels_local,&
                       recvbuf=state%parcels%numparcels_global,&
                       count=1,&
                       datatype=MPI_PARCEL_INT,&
                       op=MPI_SUM,&
                       comm=state%parallel%monc_communicator,&
                       ierror=ierr)

    !print *, "dx, dy, dz=", dx, dy, dz

    !start and end indices of the bit of grid belonging to that process
    xstart=state%local_grid%halo_size(3)+1
    ystart=state%local_grid%halo_size(2)+1
    zstart=state%local_grid%halo_size(1)+1

    xstop=nx-state%local_grid%halo_size(3)!-1
    ystop=ny-state%local_grid%halo_size(2)!-1
    zstop=nz-state%local_grid%halo_size(1)-1

    !call MPI_Barrier(state%parallel%monc_communicator,ierr)

    !print*, "rank, xstart, xstop=", state%parallel%my_rank,x_coords(xstart), x_coords(xstop+1)
    !print*, "rank, ystart, ystop=", state%parallel%my_rank,y_coords(ystart), y_coords(ystop+1)
    !print*, "rank, zstart, zstop=", state%parallel%my_rank,z_coords(zstart), z_coords(zstop+1)

    !call MPI_Barrier(state%parallel%monc_communicator,ierr)

    !call MPI_Finalize(ierr)

    !error stop



    !loop over every cell and put n_per_cell parcels in it
    n=0
    do i=xstart,xstop
      dx = x_coords(i+1)-x_coords(i)
      ddx=dx/n_per_dir
      do j=ystart,ystop
        dy = y_coords(j+1)-y_coords(j)
        ddy=dy/n_per_dir
        do k=zstart,zstop
          dz = z_coords(k+1)-z_coords(k)
          ddz=dz/n_per_dir

          do ii=1,n_per_dir
            do jj=1,n_per_dir
              do kk=1,n_per_dir
                n=n+1
                xs=x_coords(i)+ddx/2+(ii-1)*ddx
                ys=y_coords(j)+ddy/2+(jj-1)*ddy
                zs=z_coords(k)+ddz/2+(kk-1)*ddz

                state%parcels%x(n) = xs
                state%parcels%y(n) = ys
                state%parcels%z(n) = zs
                state%parcels%vol(n) = 1.0
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo



    ! do n=1,8
    !   print *, state%parcels%x(n), state%parcels%y(n), state%parcels%z(n)
    ! enddo

    !print*, maxval(state%parcels%x(1:n)), xmax
    !print*, maxval(state%parcels%y(1:n)), ymax
    !print*, maxval(state%parcels%z(1:n)), zmax


    !print *, n, nparcels
    if (n .ne. nparcels) error stop "incorrect parcel numbers"

    !set up q values
    do n=1,state%parcels%numparcels_local
      do q=1,state%parcels%qnum
        state%parcels%qvalues(q,n) = q
      enddo
    enddo

    if (state%parallel%my_rank .eq. 0) print*, "parcels initialised"

    call timer_stop(handle)

  end subroutine

end module
